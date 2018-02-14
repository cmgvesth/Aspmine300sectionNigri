"""
Module for pre work with the database.

mysqlSmChecker will check for missing data regarding SM proteins in gff, proteins, protein_has_ipr and smurf tables.

tmpSmBiTable creates a temporary table of joined blast and smurf for the subsequent createBidirSmurf function/

createBidirSmurf finally creates the bidirectional blast table for smurf entries.


Todo:
* Command line usage can be enable as e.g. python3 createSmBiblastTable blblast_base smurf_bidir_hits_for_set.
* Needs argparse for two main components:
* 1. Checking data
* 2. running tables
"""


import sys
import logging
import bioSlim3 as bio
import argparse
import MySQLdb as mdb

import misc

with open("config.txt") as c:
	config = misc.readConfig(c.readlines())



#################################
# Functions to check mysql data integrity

def parseOrgProt(l):
	"""
	Creates a dictionary of org_ids (key) and protein_ids (list of values). Handy for our table structure to quickly get all protein ids for each org.
	"""
	d = {}
	for org, prot in l:
		if str(org) not in d:
			d[str(org)] = []
		d[str(org)].append(str(prot))
	return(d)


def mysqlSmChecker(orgSet, logfile):
	"""
	Checking tables made easy!
	Please provide an orgSet and a logfile name to check the gff, proteins, protein_has_ipr and smurf tables for matching protein_ids of your selected organisms.
	"""

	logging.basicConfig(filename=logfile,level=logging.DEBUG)

	# Getting org_ids from organism table
	handle = bio.dbFetch("""SELECT org_id, name, real_name, section  FROM organism WHERE name IN ('%s')""" % "','".join(orgSet))

	orgIds = {}
	for org_id, name, real_name, section in handle:
		try:
			int(org_id)
			orgIds[name] = str(org_id)
		except Exception as e:
			print("Species %s does not have a valid org_id" % name)
			logging.error("Organism %s does not have a valid org_id" % name)

	set(orgIds.keys())
	org_ids = set(orgIds.values())




	# Getting data from server
	print("Downloading smurf data")
	smurfRaw = bio.dbFetch("""SELECT org_id, sm_protein_id, sm_short, clust_backbone
	FROM smurf
	WHERE org_id IN ('%s')
	GROUP BY org_id, clust_backbone, sm_protein_id""" % "','".join(org_ids) )

	smurf = parseOrgProt([(org, prot) for org, prot, sm, clust_backbone in smurfRaw])
	smurfSMonly = parseOrgProt([(org, prot) for org, prot, sm, clust_backbone in smurfRaw if sm != "none"])

	# Special part to check sm_proteins per backbone

	smurf_bb = {}
	for org, prot, sm, clust_backbone in smurfRaw:
		if str(org) not in smurf_bb:
			smurf_bb[str(org)] = {}
		if str(prot) not in smurf_bb[str(org)]:
			smurf_bb[str(org)][str(prot)] = []
		smurf_bb[str(org)][str(prot)].append(str(clust_backbone))

	print("Downloading proteins data")
	proteinsRaw = bio.dbFetch("""SELECT org_id, prot_seqkey FROM proteins WHERE org_id IN ('%s') GROUP BY org_id, prot_seqkey""" % "','".join(org_ids) )

	proteins = parseOrgProt(proteinsRaw)

	print("Downloading protein has ipr data")
	protein_has_iprRaw = bio.dbFetch("""SELECT org_id, protein_id FROM protein_has_ipr WHERE org_id IN ('%s') GROUP BY org_id, protein_id""" % "','".join(org_ids) )

	protein_has_ipr = parseOrgProt(protein_has_iprRaw)

	# gff
	print("Downloading gff data")
	gff_raw = bio.dbFetch("""SELECT org_id, gff_protein_id FROM gff WHERE org_id IN ('%s')""" % "','".join(org_ids) )

	gff = parseOrgProt(gff_raw)

	# Checking data for consistency

	print("Starting data check, see %s for further info" %logfile)

	for name, org_id in orgIds.items():
		org_id = str(org_id)

		try:
			smurf_bb[str(org_id)]
			try:
				smurf_bb[str(org_id)].values()
				multiBbPerProt = [item for item in smurf_bb[str(org_id)].values() if len(item) >1]
				if multiBbPerProt:
					multiBbPerProt = list(set([item for sublist in multiBbPerProt for item in sublist]))
					logging.warning("%s Identical sm_protein_ids for the following cluster backbones: %s" % (name,",".join(multiBbPerProt)))
					print("%s Identical sm_protein_ids for the following cluster backbones: %s" % (name,",".join(multiBbPerProt)))

			except Exception as e:
				logging.info("Smurf data ok for %s" % name)

		except Exception as e:
			logging.error("No smurf files for %s", name)

		try:
			smurf[org_id]
			protIdCheck = set(smurf[org_id]) - set(proteins[org_id])
			try:
				gffCheck = set(smurf[org_id]) - set(gff[org_id])
				if gffCheck:
					logging.warning("Missing gff entries for %s", name)
			except Exception as e:
				logging.error("Could not get gff files for %s", name)
			if protIdCheck:
				logging.warning("Some SM ids have not been found in proteins")
				logging.warning(protIdCheck)
			else:
				logging.info("%s has correct protein_ids in smurf", name)

			iprSMCheck= set(smurfSMonly[org_id]) - set(protein_has_ipr[org_id])
			if iprSMCheck:
				logging.error("Interpro entries are missing for major SM proteins")
			else:
				iprAllCheck= set(smurf[org_id]) - set(protein_has_ipr[org_id])
				logging.info("%s: %s out of %s SM proteins do not have annotations" % (name, len(iprAllCheck), len(set(smurf[org_id]))) )

		except Exception as e:
			logging.warning("Organism: %s, org_id: %s does not contain smurf entries" %(name, org_id))


##################################
# FUNCTIONS ON MYSQL TABLE CREATION

def tmpSmBiTable(smtable = "null"):
	"""
	This function creates a smtable of biblast hits between gene clusters.
	A biblast smtable must exist for the secondary metabolite analysis pipeline
	in order to run.
	"""

	if smtable == "null":
		raise ValueError("smtable argument needs to be provided. Try e.g. tmpSmBiTable('I_Flavi_biblast')")

	db = mdb.connect(host=config['host'], user=config['user'], passwd=config['passwd'], db=config['db'])
	cursor = db.cursor()

	# Checking if tables are available:
	query = "SHOW TABLES"
	cursor.execute(query)
	tables = [item[0] for item in  cursor.fetchall()] # reminder to self: If I fetch from mysql like this there will always be a tuple returned


	if smtable not in tables:
		db.close()
		raise ValueError("Requested smtable %s is not available in database" % smtable)

	if "smurfTemp" in tables:
		print("Table smurfTemp already exists. If you performed a run on the same set before, keep it.\
		If it is on another set or you don't know, delete it.\n")

		smurfTempAns = input("Delete smurfTemp table and generate new one? (y/n)\n").lower()
		if smurfTempAns.startswith('y'):
			print("Deleting old temporary table 1 (smurfTemp)")
			query = """DROP TABLE smurfTemp;"""
			cursor.execute(query)
			print("Done")

			print("Creating new smurfTemp table")
			query = """
			CREATE TABLE smurfTemp AS
			SELECT smurf.*, organism.name
			FROM smurf LEFT JOIN organism USING(org_id);
			CREATE INDEX i_name ON smurfTemp (name);
			CREATE INDEX i_protein ON smurfTemp (sm_protein_id);"""

			cursor.execute(query)
			print("Done\n")
		else:
			print("Keeping old smurfTemp table\n")

	if "smurf_bidir_hits_tmp1" in tables:
		print("Table smurf_bidir_hits_tmp1 already exists. If you performed a run on the same set before, keep it. If it is on another set or you don't know, delete it.")
		smurf_bidir_hits_tmp1_Ans = input("Delete smurf_bidir_hits_tmp1 table and generate new one? (y/n)\n").lower()
		if smurf_bidir_hits_tmp1_Ans.startswith('y'):
			print("Deleting temporary smtable smurf_bidir_hits_tmp1\n")
			query = """DROP TABLE smurf_bidir_hits_tmp1;"""
			cursor.execute(query)
			print("Done")

			query = """
			CREATE TABLE smurf_bidir_hits_tmp1 AS
			  SELECT smurfQ.org_id AS q_org, smurfQ.sm_protein_id AS q_protein_id, CONCAT(smurfQ.org_id, '_' , smurfQ.clust_backbone,'_', smurfQ.clust_size) AS q_clust_id, bidir.h_org, bidir.h_seqkey, bidir.pident, bidir.q_cov, bidir.h_cov
			  FROM %s AS bidir
			  JOIN smurfTemp AS smurfQ ON bidir.q_org = smurfQ.name AND bidir.q_seqkey = smurfQ.sm_protein_id;""" % smtable

			print("Creating new temporary smtable smurf_bidir_hits_tmp1")

			print("Executing query")
			print(query)
			print("Yes, this will take a while....")
			try:
				cursor.execute(query)
			except Exception as e:
				raise print("Could not run smurf bidir hits query for temporary smtable, check your connection to the database")

			indexQueries = ["CREATE INDEX i_h_org ON smurf_bidir_hits_tmp1 (h_org);",
			"CREATE INDEX i_h_protein ON smurf_bidir_hits_tmp1 (h_seqkey);"]

			print("Creating indexes")
			for query in indexQueries:
				cursor.execute(query)

			print("Done")

			print("Done")
		else:
			print("Keeping smurf_bidir_hits_tmp1")
		db.close()
	# biblast_ID50_SC130

def bidirExec(cursor,smtable):

	query = """CREATE TABLE %s AS SELECT bidir.q_org, bidir.q_protein_id, bidir.q_clust_id, smurfH.org_id AS h_org, smurfH.sm_protein_id AS h_protein_id, CONCAT(smurfH.org_id, '_' , smurfH.clust_backbone,'_', smurfH.clust_size) AS h_clust_id, bidir.pident, bidir.q_cov, bidir.h_cov
	FROM smurf_bidir_hits_tmp1 AS bidir JOIN smurfTemp AS smurfH ON bidir.h_org = smurfH.name AND bidir.h_seqkey = smurfH.sm_protein_id;""" % smtable

	print("Executing query")
	print(query)

	cursor.execute(query)
	print("Done")


	indexQueries = ["CREATE INDEX i_qorg ON %s (q_org);" % smtable,
	"CREATE INDEX i_horg ON %s (h_org);" % smtable,
	"CREATE INDEX i_qprotein ON %s (q_protein_id);"% smtable,
	"CREATE INDEX i_q_clust_id ON %s (q_clust_id);"% smtable,
	"CREATE INDEX i_h_clust_id ON %s (h_clust_id);"% smtable]

	print("Creating indexes")
	for query in indexQueries:
		cursor.execute(query)
	print("Done")


def createBidirSmurf(smtable = "none"):
	print("Starting function for final smurf_bidir_hits\n")
	db = mdb.connect(host=config['host'], user=config['user'], passwd=config['passwd'], db=config['db'])
	cursor = db.cursor()

	query = "SHOW TABLES"
	cursor.execute(query)
	tables = [item[0] for item in  cursor.fetchall()]

	if smtable in tables:
		print("Table %s already exists. If you performed a run on the same set before, keep it.\
		If it is on another set or you don't know, delete it.\n" % smtable)

		smurfTempAns = input("Delete %s table and generate new one? (y/n)\n" % smtable).lower()
		if smurfTempAns.startswith('y'):
			print("Deleting temporary smtable %s" % smtable)
			query = """DROP TABLE %s;""" % smtable
			cursor.execute(query)

			bidirExec(cursor, smtable)
		else:
			print("Keeping %s \n" % smtable)
	else:
		bidirExec(cursor, smtable)

	print("Done")
	db.close()
