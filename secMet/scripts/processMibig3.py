"""
Module to dereplicate our SM proteins with the mibig data.

Todo:
* Make it a standalone script with main?
"""

import os
import csv
import sys
import urllib.request
import tarfile
# sys.path.append("..")
import bioSlim3 as bio
import argparse
from Bio import SeqIO
from collections import Counter
import pandas as pd

def dlMibig():
	"""
	Creates folder to download and extract mibig files. From https://mibig.secondarymetabolites.org/mibig_gbk_1.3.tar.gz.
	"""
	os.getcwd()
	print("Creating directory for mibig files")
	if "mibigGbkFiles" not in os.listdir():
		os.mkdir("mibigGbkFiles")

	os.chdir("mibigGbkFiles")
	try:
		print("Downloading")
		urllib.request.urlretrieve("https://mibig.secondarymetabolites.org/mibig_gbk_1.3.tar.gz", "mibig.tar.gz")
		print("Extracting")
		try:
			tar = tarfile.open("mibig.tar.gz")
			tar.extractall()
			tar.close()
			os.chdir("..")
		except Exception as e:
			print("Cannot unzip gbk files")

	except Exception as e:
		print("Cannot download database")




def dlSmurfProteins(orgSet):
	""" Download backbones of smurf database for annotation by mibig
	needs orgSet for join with
	Returns tuples of (name,seq)\n
	orgSet must be a character vector of jgi names"""
	query = """
	SELECT smurf.org_id, smurf.sm_protein_id, proteins.prot_seq FROM smurf
	JOIN organism ON organism.name IN ('%s') AND smurf.org_id = organism.org_id
	JOIN proteins ON smurf.sm_short != 'none' AND smurf.org_id = proteins.org_id
	AND smurf.sm_protein_id = proteins.prot_seqkey;
	""" % "','".join(orgSet)

	print("Query to download smurf entries")
	print(query)

	sProts = bio.dbFetch(query)

	sp = [(str(org) + '_' + str(prot), seq.decode("UTF-8")) for org, prot, seq in sProts]

	print(sp[0:10])

	return(sp)



def processMibig(mibigDir, selOrgs = ["Aspergillus", "Penicillium"]):
	"""Provide the mibig db as genbank file
	returns SeqIO records of gene clusters
	ONly looking for SMGC from Aspergillus and Penicillium"""

	# Filtering for gbk files from mibig
	files = [file for file in os.listdir(mibigDir) if file.startswith("BGC")]

	geneClusters = []
	failedOnes = []


	for file in files:
		try:
			with open(mibigDir + "/" +file) as tempFile:
				clusters = [rec for rec in SeqIO.parse(tempFile, 'genbank')]
				if selOrgs != "all":
				# Filtering gene clusters for Aspergillus organism
					for cluster in clusters:
						if any(org in cluster.annotations["organism"] for org in selOrgs):
							geneClusters.append(cluster)
						else:
							continue
				else:
					for cluster in clusters:
						geneClusters.append(cluster)
		except:
			print("Failed to read" + file)
			failedOnes.append(file)

	print ("%s gene cluster entries out of %s genbank files could be loaded" % (len(geneClusters), len(files)))
	print ("%s files could not be read" % len(failedOnes))
	print (" Failed files:")
	print (failedOnes)

	return(geneClusters)


def getLargestProtein(seqRec):
	"""
	Use on seqIO object to get largest protein sequence.
	We want to get the largest protein sequence since SM proteins like NRPS, PKS, etc. are the largest ones in secondary metabolic gene clusters.
	"""
	bb = []
	if seqRec.features:
		for feature in seqRec.features:
			if feature.type == "CDS":
				if 'translation' in feature.qualifiers:
					bb.append(feature.qualifiers['translation'][0])
	backbone = max(bb, key = len)
	return(backbone)



#
def writeMibigFormatted(geneClusters, addPath = ''):
	"""
	Get largest proteins from each cluster and write them to mibigDf.fasta
	Sometimes the files just contain domains only and not the full protein.
	That's ok, we will still get the best hit at a later stage.

	Output:
	* mibigDb.fasta (containing mibig entries in blastable format)
	* translateBgc.txt (mibig entries with their corresponding blastable id)
	"""
	finDb = {}
	translationTable = []
	counter = 1
	for rec in geneClusters:

		comp = rec.description
		org = rec.annotations['organism']
		seq = getLargestProtein(rec)
		gcId = "gc_"+str(counter)

		finDb[gcId] = seq
		translationTable.append((comp, org ,gcId))
		counter = counter + 1

	# WRITING SEQUENCES #
	print("Writing sequences to %s" % addPath + 'mibigDb.fasta')
	bio.dictToFasta(finDb, addPath + "mibigDb.fasta")

	# WRITING TRANSLATION TABLE #
	print("Writing translation table to %s" % addPath + 'translateBgc.txt')
	with open(addPath + "translateBgc.txt", "w") as tfile:
		out = csv.writer(tfile, delimiter=';')
		for row in translationTable:
			out.writerow(row)

def processBlastResult(blastFile, translationFile, pident_cutoff = 95, q_coverage_cutoff = 90, s_coverage_cutoff = 90, aspmineIdsOnly = True):
	"""
	Processes the blast results of mibig to find the best hit in our data.
	Needs the name of the blast file and the name of the generated translation table from writeMibigFormatted.

	"""
	bheader = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send","evalue", "bitscore", "qlen", "slen"]

	try:
		blast = pd.read_csv(blastFile, header = None, names = bheader , sep = "\t")
	except Exception as e:
		raise
		print("Check your blast file")

	blast['q_cov'] =  (blast['qend'] - blast['qstart']) / blast['qlen']*100
	blast['s_cov'] =  (blast['send'] - blast['sstart']) / blast['slen']*100

	subBlast = blast[(blast.pident > pident_cutoff) & (blast.s_cov > s_coverage_cutoff) & (blast.q_cov > q_coverage_cutoff)]
	# "/home/seth/asptoolbox/secMet/mibig/translateBgc.txt"
	annotation = pd.read_csv(translationFile, names = ["compound", "org", "gcId"], sep = ";")

	subBlastAnnoteted = pd.merge(subBlast, annotation,  how='left', left_on=['qseqid'], right_on = ['gcId'])

	# Getting the best match
	inds = subBlastAnnoteted.groupby(['compound'])['pident'].idxmax()
	subBlastFinal = subBlastAnnoteted.ix[inds]

	if aspmineIdsOnly: # Need to do this so I can use the module for ncbi identifiers...
		subBlastFinal['org_id'] = [item[0] for item in subBlastFinal.sseqid.str.split("_")]
		subBlastFinal['protein_id'] = [item[1] for item in subBlastFinal.sseqid.str.split("_")]

	return(subBlastFinal)
