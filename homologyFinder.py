#!/usr/bin/python

# python homologyFinder.py -host <host> -user <user> -passwd <passwd> -btable <biblast_tablename> -htable <homology_tablename> -odir <output_dir> -spfile /path/to/species_file.txt 

# DESCRIPTION: 
# This program single-link blast hits into homologous protein families (hfams) creating a 
# homology table on Aspmine. It further couples the both Interpro and GO annotations to the proteins
# in the homology table, creating a homology_InterPro table and a homology_GO table respectively. 
# This program retrieves its data from the Aspmine database. Normally the blast files have been 
# filtered prior to this step and combined into one table called `biblast`. 
# Usual cutoff blast alignment values are: identity >= 50%, (q_cov + h_cov) >= 130% and the hits have to be reciprocal.    
# This program can either create a new `homology table` or append to an existing table based on species inputted
# in the commandline.

# INPUT:
# 1. species_file.tsv
# A file containing the JGI species names - one on each line.
# E.g.
# Aspamy1
# Aspcal1
# Aspfu1

# OUTPUT:
# 1. homology table
# A table with the single-linked protein families based on biblast records
# hfam   org_id   org_name   protein_id
# 1	     27	      Aspoch1	 100261
# 2	     27	      Aspoch1	 100324
# 3	     27	      Aspoch1	 10057
# 4	     27	      Aspoch1	 10174
# INDEXES: `i_hfam_org_prot` (`hfam`,`org_name`,`protein_id`) and `i_org_prot` (`org_name`,`protein_id`)

# 2. homology_IPR table
# A table where the InterPro annotation is appended to the proteins in the homology table
# hfam   org_id   org_name   protein_id   ipr_id   ipr_desc
# 1	     27	      Aspoch1	 100261	      NULL	   NULL
# 2	     27	      Aspoch1	 100324	      NULL	   NULL
# 3	     27	      Aspoch1	 10057	      NULL	   NULL
# 4	     27	      Aspoch1	 10174	      NULL	   NULL

# 3. homology_GO table
# A table where the GO annotation is appended to the proteins in the homology table
# hfam   org_id   org_name   protein_id   go_term_id   go_name   go_termtype
# 1	     27	      Aspoch1	 100261	      NULL	       NULL      NULL
# 2	     27	      Aspoch1	 100324	      NULL	       NULL      NULL
# 3	     27	      Aspoch1	 10057	      NULL	       NULL      NULL
# 4	     27	      Aspoch1	 10174	      NULL	       NULL      NULL

# 4. homologyFinder.log
# A log file with all the screen outputs - located in the output directory 

#--------------------------------------------------------
# IMPORTS
#--------------------------------------------------------
import os, datetime, getpass
import sys
from sys import argv

from collections import defaultdict
import datetime
from inspect import currentframe, getframeinfo
import inspect
import itertools
import time

import operator
import errno

import shutil
import subprocess

import getopt, argparse, re, glob
from argparse import ArgumentParser

import logging

import MySQLdb as mdb

#--------------------------------------------------------
# SUBFUNCTIONS
#--------------------------------------------------------
""" FORMAT ARGUMENT HELPER """
class SmartFormatter(argparse.HelpFormatter):
	width = 100
	def _split_lines(self, text, width):
		if text.startswith('R|'):
			return text[2:].splitlines()
		return argparse.HelpFormatter._split_lines(self, text, width)

""" FORMAT ARGUMENT PARSER """
""" Allows help display when no arguments are given """
class CustomArgumentParser(argparse.ArgumentParser):
	def error(self, message):
		print("#--------------------------------------------------------------")
		print("# HELP:")
		print("#--------------------------------------------------------------")
		self.print_help()
		print("#--------------------------------------------------------------")
		print("# ERROR: %s" %message)
		print("#--------------------------------------------------------------")
		sys.exit()

""" CONNECT TO DATABASE """
def connect_db(dbname, linenumber):
	try:
		db = mdb.connect(host=args.host, user=args.user, passwd=args.passwd, db=args.dbname)
	except mdb.Error as e:
		sys.exit("# ERROR line %s - %d: %s" % (linenumber, e.args[0],e.args[1]))
	try:
		cursor = db.cursor()
	except mdb.Error as e:
		sys.exit("# ERROR line %s - %d: %s" % (linenumber, e.args[0],e.args[1]))
	return db, cursor

""" COMBINE EXECUTE AND FETCH ALL """
def executeQuery(cursor, query):
	(columns, result) = ([],[])
	try:
		cursor.execute(query)
		result = cursor.fetchall()
		if result:
			if len(result) > 0:
				columns = map(lambda x:x[0], cursor.description) 	
	except mdb.Error as e:
		sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))
	return columns, result

""" FETCH ALL MEMBERS FROM HFAM (homoTable) AND DELETE EXISTING ONES """
def hfamMembers_homoTable(hfam_homoTable_collected, homoTable, args):
	hfamMembers_homoTable = list()

	format_hfam = ', '.join(str(hfam_entry) for hfam_entry in hfam_homoTable_collected)	# Convert integers into string format

	# Retrieve all members in hfam
	query_hfam_data = ("SELECT * from %s WHERE hfam IN (%s);" % (homoTable, format_hfam))
	(column_names, hfam_data) = executeQuery(cursor, query_hfam_data)
	hfamMembers_homoTable = list(hfam_data)

	# Delete members with the specific hfams
	query = "DELETE from %s WHERE hfam IN (%s);" % (homoTable, format_hfam)
	(column_names, delete_data) = executeQuery(cursor, query)
	db.commit()

	return (hfamMembers_homoTable)


""" CREATE NEW COLLECTED FAMILY (hfam) """
def createNewFamily(hfamMembers, upload_counter, new_hfam, new_family_to_values2insert):

	for member in hfamMembers:
		# Append entries to existing cluster
		values = (int(new_hfam), member[1], member[2], member[3])
		new_family_to_values2insert.append(values)
		upload_counter += 1

	return(upload_counter, new_family_to_values2insert)

#--------------------------------------------------------
# ARGUMENTS, SETUPS AND PRINT
#--------------------------------------------------------
""" TIME """
startTimet_1 = datetime.datetime.now()
startTime = datetime.datetime.now().time() # record runtime
today = datetime.date.today()
now = datetime.datetime.now()

""" DATABASE """
parser = CustomArgumentParser(formatter_class=SmartFormatter, usage='%(prog)s -dbname [database name]')
parser.add_argument("-dbname", required=False, default = "aspminedb", help="Database name")
parser.add_argument("-host", required=True, help="Host name")
parser.add_argument("-user", required=True, help="User name")
parser.add_argument("-passwd", required=True, help="Password")

""" INPUTS """
parser.add_argument("--biblasttable", "-btable", required=True, type = str, default = '', help="Input BLAST table used for the paralog and homolog analysis. Ex. biblast_ID[custom]_SC[custom]. ID = min. alignment identity, SC = min. collected alignment coverage [SUMCOV = min.(q_cov + h_cov)].")
parser.add_argument("--homotable", "-htable", required=True, type = str, default = '', help="Input homolog table containg homologous protein family, species and protein names. If present, new species homologs will be appended to the existing table. Ex. 'homoPF_SC[custom]_proteins'")
parser.add_argument("--output_dir", "-odir",  required=True, default="", type = str, help="Name and path of the output directory.") 
""" SPECIES SELECTION """
parser.add_argument("--species", "-sp", nargs = '*',  required=False, default=[], action='store', help="List of species to analyse. Please insert JGI species names. Do NOT use any comma or quotes.")
parser.add_argument("--speciesfile", "-spfile", required=False, type = str, default = "", help="The absolute path and name of the file containing one JGI species name per line.") 
parser.add_argument("-all", required=False, action='store_true', help="Create homolog table based on all species present in the biblast table")
""" FLAGS """
parser.add_argument("-nogo", required=False, action='store_true', help="Do not create homolog table including GO terms")
parser.add_argument("-noipr", required=False, action='store_true', help="Do not create homolog table including Interpro descriptions")

""" PARSE ARGUMENTS """
args = parser.parse_args()
# database
dbname = args.dbname
host = args.host
user = args.user
passwd = args.passwd
# inputs
biblastTable = args.biblasttable
homoTable = args.homotable
output_dir = args.output_dir
# species
species = args.species 
all_species = args.all
speciesfile = args.speciesfile
# flags
nogo = args.nogo
noipr = args.noipr

#------------------------------------------------------------------
# Creating output directory if it doesn't exist
#------------------------------------------------------------------
# The code has to be located here because the logging file is written to the output folder
if os.path.isdir(output_dir) == False:
	print "# WARNING: The output directory did not exists - The program will create: %s" % output_dir
	try:
		os.makedirs(output_dir)
	except OSError as exc:
		# If the folders already exist - then pass
		if exc.errno == errno.EEXIST and os.path.isdir(output_dir):
			pass
		else: 
			raise

#--------------------------------------------------------
# LOGGING
#--------------------------------------------------------
# define a Handler which writes INFO messages or higher to the sys.stderr
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-10s %(levelname)-10s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='%s/homologyFinder.log' %(output_dir),
                    filemode='w')

console = logging.StreamHandler() # define the consol
console.setLevel(logging.INFO) # show the same information in the consol
formatter = logging.Formatter('%(name)-10s: %(levelname)-10s %(message)s') # simple consol format
console.setFormatter(formatter) # tell the handler to use this format
logging.getLogger('').addHandler(console) # add the handler to the root logger


#------------------------------------------------------------------
# Check commandline arguments
#------------------------------------------------------------------

if species == [] and not all_species and speciesfile == "":
	sys.exit("# ERROR: Please select either all species (-all) or add species to input list (-sp X Y Z) or as a file (-spfile)")
if ((species != [] and all_species) or (species != [] and speciesfile != "") or (speciesfile != "" and all_species)):
	sys.exit("# ERROR: Please select only one of the flags (-all or -sp X Y Z or -spfile <file.txt>)")
if speciesfile != "":
	if not os.path.isfile(speciesfile):
		logging.error('The species files (-spfile) did not exists: %s' %speciesfile)
		sys.exit()

#------------------------------------------------------------------
# Print argument values to screen
#------------------------------------------------------------------
# Start information on consol and log file
logging.info('--------------------------------------------------------------')
logging.info('INFO:')
logging.info('--------------------------------------------------------------')
logging.info('python ' + ' '.join(argv))
logging.info('DATE: %s'%now.strftime("%a %b %d %Y %H:%M"))
logging.info('USER: ' + getpass.getuser())
logging.info('CWD: %s' %os.getcwd())
logging.info('--------------------------------------------------------------')
logging.info('ARGUMENTS:')
logging.info('--------------------------------------------------------------')
logging.info('Database: %s' %dbname)
logging.info('Homolog table: %s' %homoTable)
if species != []:
	logging.info('List of input species: %s' %species)
if all_species:
	logging.warning('List of input species: %s =  will be found in %s' %(species, homoTable))
if speciesfile != "":
	logging.info('Input species file: %s' %speciesfile)
logging.info('Output directory: %s' %output_dir)
logging.info('--------------------------------------------------------------')

#--------------------------------------------------------
# FINDING MISSING ORGANISMS IN BIBLAST AND HOMOLOG TABLES
#--------------------------------------------------------
create_table = False
homoTable_orgs = list()
missing_orgs = list()
biblastTable_orgs = list()
all_possible_orgPairs = list()
diff_biblast_allPossible = list()

""" EXTRACT BIBLAST ORGS """
logging.info('Retrieving information from table: %s' %biblastTable)
db, cursor = connect_db(args, inspect.stack()[0][2])
if cursor.execute("Show tables LIKE '%s'" %biblastTable):
	# Extract q_org and h_org names from biblast table
	biblast_qhorg_query = "SELECT DISTINCT q_org, h_org FROM %s" %biblastTable
	(column_names, biblastTable_qhorgs) = executeQuery(cursor, biblast_qhorg_query)

	for org_pair in biblastTable_qhorgs:
		if org_pair[0] not in biblastTable_orgs:
			biblastTable_orgs.append(org_pair[0])
		if org_pair[1] not in biblastTable_orgs:
			biblastTable_orgs.append(org_pair[1])
else:
	logging.error('%s did not exist. Please rerun with a new biblast table.' %biblastTable)


""" EXTRACT HOMOTABLE ORGS AND MAX(HFAM) """
if homoTable != "":
	logging.info('Retrieving information from table: %s' %homoTable)
	if cursor.execute("SHOW TABLES LIKE '%s';" % homoTable):
		# Extract org_names from homolog table
		homo_org_query = "SELECT DISTINCT org_name FROM %s" %homoTable
		(column_names, homoTable_orgs) = executeQuery(cursor, homo_org_query)
		homoTable_orgs = map(' '.join, homoTable_orgs) 	# Convert list of tuples to list of strings

		# Find the max homolog cluster number (hfam)
		(column_names, max_hfam) = executeQuery(cursor, "SELECT MAX(hfam) FROM %s;" %homoTable)
		if max_hfam[0][0] == None:
			new_hfam = int(1)
		else:
			new_hfam = int(max_hfam[0][0]+1)
	else:
		logging.info('%s did not exist - creating new table' %homoTable)
		create_table = True
		new_hfam = int(1)
else:
	create_table = True
	new_hfam = int(1)


""" CHECK COMMANDLINE INPUT SPECIES """
if speciesfile != "":
	species = list()
	sp_file = open(speciesfile, "r")
	species = [line.strip() for line in sp_file]
	sp_file.close()

if species != []:
	# Check if all input species are in biBLAST table
	org_not_in_biblast_orgs = list(set(species) - set(biblastTable_orgs))
	if len(org_not_in_biblast_orgs) > 0:
		logging.error('Please recreate %s including species from input: %s' % (biblastTable, org_not_in_biblast_orgs))
		sys.exit()

if all_species:
	# IF selected all orgs from biblast table
	species = biblastTable_orgs

# Combine all input species with the orgs from the homolog table
# and create all possible orgPairs
input_homoTable_species = set(homoTable_orgs + species)
all_possible_orgPairs = tuple(itertools.product(set(input_homoTable_species), set(input_homoTable_species)))

# Check that all possible orgPairs are in biblast table
diff_biblast_allPossible = list(set(all_possible_orgPairs) - set(biblastTable_qhorgs))
if len(diff_biblast_allPossible) > 0:
	logging.error("These organism pairs are missing to make a complete 'all vs. all' single linkage:")
	for dif in diff_biblast_allPossible:
		logging.error('%s' % str(dif))
	sys.exit()

""" RETRIEVE ORG ID FROM organism TABLE """
# Creating lookup lists to be able to change commandline 
# input organisms to one common organism ID for further process
logging.info('Retrieving information from table: organism')
try:
	# Extract organism ID, name and real name from database
	organism_query = "SELECT org_id, name, real_name, section FROM organism"
	cursor.execute(organism_query)
	organism_name = cursor.fetchall()
except mdb.Error, e:
	logging.error("%d: %s" % (e.args[0],e.args[1]))
	db.close
	sys.exit()
db.close

# Restructure organism data from database
orgname_to_id = dict()
for element in organism_name:
	orgname_to_id[element[1]] = element[0]


""" FIND ORGS TO BE APPENDED TO HOMOTABLE """
# Find missing q_orgs from homolog table
if homoTable_orgs:
	missing_orgs = list(set(species) - set(homoTable_orgs))
else:
	missing_orgs = species

logging.info('Species to be appended in table %s: %s' %(homoTable, missing_orgs))

#--------------------------------------------------------
# CREATE NEW MYSQL HOMOLOG TABLE FOR SINGLE LINKAGE
#--------------------------------------------------------

if create_table:
	logging.info('Creating table: %s' % homoTable)
	db, cursor = connect_db(args, inspect.stack()[0][2])
	query = ("""CREATE TABLE %s (
		`hfam` int(100) NOT NULL,
		`org_id` int(100) NOT NULL,
		`org_name` varchar(100) NOT NULL,
		`protein_id` int(100) NOT NULL,
		unique `i_hfam_org_prot` (`hfam`, `org_name`, `protein_id`),
		KEY `i_org_prot` (`org_name`, `protein_id`)
		) ENGINE=MyISAM DEFAULT CHARSET=latin1""") % homoTable
	
	executeQuery(cursor, query)
	db.commit()

	if not cursor.execute("SHOW TABLES LIKE '%s';" % homoTable):
		logging.error("%s was not created" % homoTable)
		db.close()
		sys.exit()
	db.close()

#--------------------------------------------------------
# APPEND, MERGE OR CREATE NEW CLUSTERS/FAMILIES
#--------------------------------------------------------

# Limit the biBLAST search to homoTable q_orgs + runned q_orgs
searchOrgs = homoTable_orgs
startTimet_2 = ""
startTimet_3 = datetime.datetime.now()

logging.info('Creating homologous protein families')
db, cursor = connect_db(args, inspect.stack()[0][2])
org_count = 0
for q_org in missing_orgs:
	if org_count != 0:
		logging.info('Iteration time: %s' %str(datetime.datetime.now()-startTimet_2))

	org_count += 1

	startTimet_2 = datetime.datetime.now()
	logging.info("----------------------------------------------------")
	logging.info('Running species: %s - %s of %s' %(q_org, org_count, len(missing_orgs)))
	logging.info("----------------------------------------------------")

	
	# Set variables 
	total_seqkey_counter = 0
	upload_counter = 0
	total_upload_counter = 0
	values2insert = []
	
	# Limit the biBLAST search to homoTable q_orgs + runned q_orgs
	searchOrgs.append(q_org)

	""" RETRIEVE ALL q_seqkeys """
	# Retrieve q_seqkey from missing q_org
	q_org_seqkey_query = "SELECT DISTINCT q_seqkey FROM %s WHERE q_org = '%s';" %(biblastTable, q_org)
	(column_names, q_org_seqkey_list) = executeQuery(cursor, q_org_seqkey_query)
	q_org_seqkey_list = list(sum(q_org_seqkey_list, ())) # list of tuples to flat list


	# Find all q_seqkey biBLAST hits (h_seqkeys) in both the input and homoTable q_orgs 
	for q_seqkey in q_org_seqkey_list:
		
		# Set and reset variables 
		total_seqkey_counter += 1
		blast_homoTable_output = list()
		hfam_homoTable_collected = list()
		hfam_values2insert_collected = list()
		collect_NULLhfam = list()

		""" RETRIEVE ALL HOMOLOGS TO q_seqkey """ 
		# Combine homoTable and biblast - Output: hfam;h_org;h_seqkey
		query_blast_homoTable = """SELECT hfam, tb.*
			FROM (
			# 2.1. Select only one column pair (here h_org, h_seqkey)
			SELECT h_org, h_seqkey
			FROM (
				# 1. Get all biblast hits to q_org/q_seqkey from db 
				SELECT q_org, q_seqkey, h_org, h_seqkey
				FROM %s
				WHERE ((q_org = '%s' AND q_seqkey = %s) OR (h_org = '%s' AND h_seqkey = %s))) ta
			# 2.2. where both q_org and h_org is in the selected orgs
			WHERE (ta.q_org IN ('%s') 
			AND ta.h_org IN ('%s'))
			) tb
			# 3.1 Join hfam from homotable where h_org/h_seqkey exists
			LEFT JOIN %s
			ON h_org = org_name AND h_seqkey = protein_id
			# 3.2 group to reduce duplicates and retrieve all potential hfams per h_org/h_seqkey
			GROUP BY hfam, h_org, h_seqkey
			;""" %(biblastTable, q_org, q_seqkey, q_org, q_seqkey, "', '".join(searchOrgs), 
				"', '".join(searchOrgs), homoTable)

		(column_names, blast_homoTable_data) = executeQuery(cursor, query_blast_homoTable) 


		# Converto output from tuples of tuples to a list of lists
		# or create empty list
		if len(blast_homoTable_data) == 0:
			blast_homoTable_output = list()
		else:
			blast_homoTable_output = blast_homoTable_data

		""" FETCH ALL HFAMS FROM homoTable """
		for blasthit in blast_homoTable_output:
			hfam_homoTable = blasthit[0]
			h_org = blasthit[1]
			h_seqkey = blasthit[2]

			# Retrieve all homoTable HFAMS
			if hfam_homoTable != None: # can hfam be collected from homologs in the query?
				if hfam_homoTable not in hfam_homoTable_collected:
					hfam_homoTable_collected.append(hfam_homoTable)
			else:
				# Collect all entries with no HFAM - Exit if not a paralog
				collect_NULLhfam.append(blasthit)


		""" FETCH ALL MEMBERS FROM HFAM AND DELETE EXISTING ONES AND CREATE NEW HFAMS"""
		members_homoTable = list()
		new_family_members = list()
		new_family_to_values2insert = list()
		collect_NULLhfam_orgID = list()
		upload_counter = len(values2insert)

		# Fetch hfam members, delete existing ones and create new uploads
		if len(hfam_homoTable_collected) > 0:
			members_homoTable = hfamMembers_homoTable(hfam_homoTable_collected, homoTable, args)

		# Include org_id to collect_NULLhfam protein members
		if len(collect_NULLhfam) > 0:
			for row_entry in collect_NULLhfam:
				collect_NULLhfam_orgID.append((row_entry[0], int(orgname_to_id[row_entry[1]]), row_entry[1], int(row_entry[2])))
		
		# Combine homoTable members with biblast hfams and hits without hfams
		new_family_members = list(set(members_homoTable + collect_NULLhfam_orgID))

		# Create new family with new hfam
		if len(new_family_members) > 0:
			(upload_counter, new_family_to_values2insert) = createNewFamily(new_family_members, upload_counter, new_hfam, new_family_to_values2insert)	

		# Remove duplicates in new_family_to_values2insert and extend values2insert
		# OBS: This might be redundant - IT IS NOT (tested)
		new_family_to_values2insert_reduced = list(set(new_family_to_values2insert))
		
		# Update values2insert and hfam count
		values2insert.extend(new_family_to_values2insert_reduced)
		new_hfam += 1
		
		""" UPLOAD TO SERVER """
		# Uploading to server
		if upload_counter >= 5000 or total_seqkey_counter == len(q_org_seqkey_list):

			total_upload_counter = total_upload_counter + upload_counter
			# Prints record number that will be inserted
			if upload_counter >= 5000 : 
				logging.info('Inserting record number %s' % total_upload_counter)
			elif total_seqkey_counter == len(q_org_seqkey_list): 
				logging.info('Inserting record number %s' % total_upload_counter)
				total_upload_counter = 0

			# Upload into table
			try:
				query =  "INSERT IGNORE INTO %s (hfam, org_id, org_name, protein_id) values(%s);" % (homoTable, ("%s," * len(values2insert[0])).rstrip(","))
				
				cursor.executemany(query, values2insert)
				# Add changes to database
			 	db.commit()	

				upload_counter = 0 	# restart counter
				values2insert = []	# Empty list of values
				new_family_to_values2insert = []
			except mdb.Error, e:
				print values2insert
				logging.error('%s load %s %d: %s' % (homoTable, q_org, e.args[0],e.args[1]))
				db.close()
				sys.exit()

""" CLOSE DATABASE """
db.close()
new_hfam += 1

if len(missing_orgs)>0:
	logging.info('Iteration time: %s' %str(datetime.datetime.now()-startTimet_2))

logging.info("----------------------------------------------------")
logging.info('Total linking time: %s' %str(datetime.datetime.now()-startTimet_3))
logging.info("----------------------------------------------------")

startTimet_4 = datetime.datetime.now()

#--------------------------------------------------------
# CREATE DUPLICATE TABLE
#--------------------------------------------------------
dupl_table = homoTable+"_dupl"

""" DROP TABLE IF EXISTS """
db, cursor = connect_db(args, inspect.stack()[0][2])
if cursor.execute("SHOW TABLES LIKE '%s';" % dupl_table):
	logging.warning('Deleting existing table: %s' %dupl_table)
	executeQuery(cursor, "DROP TABLE %s;" % dupl_table)
db.commit()

logging.info('Creating protein duplication table: %s' %dupl_table)
query = """CREATE TABLE %s (
	`hfam` int(100) NOT NULL,
	`name` varchar(100) NOT NULL,
	KEY `i_hfam_name` (`hfam`,`name`),
	KEY `name` (`name`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;""" %dupl_table

executeQuery(cursor, query) 

if not cursor.execute("SHOW TABLES LIKE '%s';" % dupl_table):
	logging.error("%s was not created" % dupl_table)
	db.close()
	sys.exit()

PROTdupl_table_query = """INSERT %s
	SELECT tb.hfam, CONCAT(ta.org_name, ":",ta.protein_id) name
	FROM (
	SELECT org_name, protein_id
	FROM %s
	GROUP BY org_name, protein_id
	HAVING COUNT(DISTINCT hfam) > 1) ta
	JOIN %s tb
	ON (ta.org_name = tb.org_name AND ta.protein_id = tb.protein_id)
	;""" %(dupl_table, homoTable, homoTable)

executeQuery(cursor, PROTdupl_table_query) 
db.close()


#--------------------------------------------------------
# SINGLE-LINK ALL HFAMS WITH DUPLICATED 
# org_name/protein_id pair IN DUPLICATE TABLE
#--------------------------------------------------------
""" CHECK ROW COUNT """
db, cursor = connect_db(args, inspect.stack()[0][2])

(column_names, row_count) = executeQuery(cursor, "SELECT COUNT(*) FROM %s;" % dupl_table)
logging.info('There are %s duplications present in %s' %(row_count[0][0], homoTable))

logging.info('Mergin hfams with common org/protein pairs')
while row_count[0][0] != 0:
	# Initiate variables
	protein_list = list()
	hfam_list = list()
	protein_list_updated = list()
	hfam_list_updated = list()
	protein_missing = list()
	hfam_missing = list()

	# Get first duplicated protein
	(column_names, dupl_protein) = executeQuery(cursor, "SELECT name FROM %s limit 1;" % dupl_table)

	# Retrieve all hfams and org/prot pairs associated with the dupl proteins hfams
	prot_hfam_query= """SELECT * FROM %s WHERE hfam IN (SELECT hfam FROM %s WHERE name = '%s')
		GROUP BY hfam, name;""" %(dupl_table, dupl_table, dupl_protein[0][0])
	(column_names, prot_hfam) = executeQuery(cursor, prot_hfam_query)

	# Save newly found org/protein pairs and their hfams to lists
	for entry in prot_hfam:
		if entry[1] not in protein_list_updated:
			protein_list_updated.append(entry[1])
		if entry[0] not in hfam_list_updated:
			hfam_list_updated.append(entry[0])
	
	# While a new protein is added find its hfams and the duplicated proteins that are in the hfams
	while set(protein_list) != set(protein_list_updated):
		missing_prots = list(set(protein_list_updated)-set(protein_list))
		protein_list = protein_list_updated
		hfam_list = hfam_list_updated

		for element in missing_prots:
			prot_hfam_query= """SELECT * FROM %s WHERE hfam IN (SELECT hfam FROM %s WHERE name = '%s')
			GROUP BY hfam, name;""" %(dupl_table, dupl_table, element)
			(column_names, prot_hfam) = executeQuery(cursor, prot_hfam_query)

			for entry in prot_hfam:
				if entry[1] not in protein_missing:
					protein_missing.append(entry[1])
				if entry[0] not in hfam_missing:
					hfam_missing.append(entry[0])

			protein_list_updated = set.union(set(protein_list_updated), set(protein_missing))
			hfam_list_updated = set.union(set(hfam_list_updated), set(hfam_missing))
	
	protein_list = set.union(set(protein_list), set(protein_list_updated))
	hfam_list =  set.union(set(hfam_list), set(hfam_list_updated))

	#--------------------------------------------------------
	# DUPLICATE DELETION
	#--------------------------------------------------------
	""" Deletion of proteins in dupl_table """
	delete_prot_query = "DELETE FROM %s WHERE name IN ('%s');" %(dupl_table, "', '".join(protein_list))
	executeQuery(cursor, delete_prot_query)
	db.commit()

	# Retrieve all members hfams in homoTable 
	ALLprotsInHfams_query = "SELECT * FROM %s WHERE hfam in (%s) GROUP BY org_name, protein_id;" %(homoTable, ", ".join([str(i) for i in hfam_list]))
	(column_names, ALLprotsInHfams) = executeQuery(cursor, ALLprotsInHfams_query)	


	""" Deletion of hfams in homoTable """	
	delete_prot_query = "DELETE FROM %s WHERE hfam IN (%s);" %(homoTable, ", ".join([str(i) for i in hfam_list]))
	executeQuery(cursor, delete_prot_query)
	db.commit()
	
	#--------------------------------------------------------
	# CREATION OF NEW HFAMS
	#--------------------------------------------------------
	new_family_to_values2insert = list()
	
	for prot_member in list(ALLprotsInHfams):
		# Append entries to existing cluster
		values = (int(new_hfam), int(prot_member[1]), prot_member[2], int(prot_member[3]))
		new_family_to_values2insert.append(values)


	# Upload into table
	try:
		query =  "INSERT IGNORE INTO %s (hfam, org_id, org_name, protein_id) values(%s);" % (homoTable, ("%s," * len(new_family_to_values2insert[0])).rstrip(","))
		
		cursor.executemany(query, new_family_to_values2insert)
		# Add changes to database
	 	db.commit()	

		new_family_to_values2insert = []	# Empty list of values
	except mdb.Error, e:
		print new_family_to_values2insert
		logging.error('%s load %s %d: %s' % (homoTable, q_org, e.args[0],e.args[1]))
		db.close()
		sys.exit()

	#--------------------------------------------------------
	# CHECKING ROW COUNTS
	#--------------------------------------------------------
	(column_names, row_count) = executeQuery(cursor, "SELECT COUNT(*) FROM %s;" % dupl_table) 
	new_hfam += 1

#--------------------------------------------------------
# DELETE EMPTY DUPLICATION TABLE
#--------------------------------------------------------
logging.info('Deleting duplication table: %s' %dupl_table)
executeQuery(cursor, "DROP TABLE IF EXISTS %s;" % dupl_table)
db.commit()
db.close()


logging.info("----------------------------------------------------")
logging.info('Deletion duplication time:%s' %str(datetime.datetime.now()-startTimet_4))
logging.info("----------------------------------------------------")


if not noipr:
	startTimet_ipr = datetime.datetime.now()
	ipr_table = homoTable+"_IPR"
	logging.info("----------------------------------------------------")
	logging.info('Creating Interpro table: %s' %ipr_table)
	logging.info("----------------------------------------------------")
	logging.info("This may take some time")

	""" DROP TABLE IF EXISTS """
	db, cursor = connect_db(args, inspect.stack()[0][2])
	if cursor.execute("SHOW TABLES LIKE '%s';" % ipr_table):
		logging.warning('Deleting existing table: %s' %ipr_table)
		executeQuery(cursor, "DROP TABLE %s;" % ipr_table)
	db.commit()

	logging.info('Creating Interpro table: %s' %ipr_table)
	
	""" CREATE TABLE """
	ipr_query= """CREATE TABLE %s
		SELECT hfam, org_id, org_name, protein_id, t1.ipr_id, ipr_desc
		FROM
		(SELECT hfam, hfam.org_id, hfam.org_name, hfam.protein_id, ipr_id 
		FROM %s as hfam
		LEFT JOIN protein_has_ipr as pip
		ON (hfam.org_id = pip.org_id AND hfam.protein_id = pip.protein_id)
		GROUP BY hfam, org_name, hfam.protein_id, ipr_id) t1
		LEFT JOIN ipr
		ON (t1.ipr_id = ipr.ipr_id)
		GROUP BY hfam, org_name, protein_id, t1.ipr_id;""" %(ipr_table,homoTable)

	executeQuery(cursor, ipr_query) 

	""" INDEX """
	if not cursor.execute("SHOW TABLES LIKE '%s';" % ipr_table):
		logging.error("%s was not created" % ipr_table)
		db.close()
		sys.exit()
	else:
		logging.info('Creating indexes - this may take a while')
		executeQuery(cursor, "CREATE INDEX `i_all` ON %s (hfam, org_name, protein_id, ipr_id);" %ipr_table)
		executeQuery(cursor, "CREATE INDEX `i_org_prot` ON %s (org_name, protein_id);" %ipr_table)
		db.commit()

	# Retrieving all species from IPR table
	logging.info('Checking that all species have InterPro annotations')
	IPRorgs_query = "SELECT DISTINCT org_name FROM %s" %ipr_table
	(column_names, IPRtable_orgs) = executeQuery(cursor, IPRorgs_query)
	db.close()

	# Create warning if species are missing InterPro annotation
	IPRtable_orgs = map(' '.join, IPRtable_orgs)
	missing_orgs = list()
	missing_orgs = list(set(homoTable_orgs)-set(IPRtable_orgs))
	if len(missing_orgs) > 0:
		logging.warning("These species do not have InterPro annotation:")
		for morg in missing_orgs:
			logging.warning("%s" %morg)

	logging.info('Finished %s - runtime :%s' %(ipr_table, str(datetime.datetime.now()-startTimet_ipr)))



if not nogo:
	startTimet_go = datetime.datetime.now()
	go_table = homoTable+"_GO"
	logging.info("----------------------------------------------------")
	logging.info('Creating GO table: %s' %go_table)
	logging.info("----------------------------------------------------")
	logging.info("This may take some time")

	""" DROP TABLE IF EXISTS """
	db, cursor = connect_db(args, inspect.stack()[0][2])
	if cursor.execute("SHOW TABLES LIKE '%s';" % go_table):
		logging.warning('Deleting existing table: %s' %go_table)
		executeQuery(cursor, "DROP TABLE %s;" % go_table)
	db.commit()

	logging.info('Creating GO table: %s' %go_table)
	
	""" CREATE TABLE """
	go_query= """CREATE TABLE %s
		SELECT hfam, org_id, org_name, protein_id, t1.go_term_id, go_name, go_termtype
		FROM
		(SELECT hfam, hfam.org_id, hfam.org_name, hfam.protein_id, go_term_id 
		FROM %s as hfam
		LEFT JOIN protein_has_go as pgo
		ON (hfam.org_id = pgo.org_id AND hfam.protein_id = pgo.protein_id)
		GROUP BY hfam, org_name, hfam.protein_id, go_term_id) t1
		LEFT JOIN go
		ON (t1.go_term_id = go.go_term_id)
		GROUP BY hfam, org_name, protein_id, t1.go_term_id;""" %(go_table, homoTable)

	executeQuery(cursor, go_query) 

	""" INDEX """
	if not cursor.execute("SHOW TABLES LIKE '%s';" % go_table):
		logging.error("%s was not created" % go_table)
		db.close()
		sys.exit()
	else:
		logging.info('Creating indexes - this may take a while')
		executeQuery(cursor, "CREATE INDEX `i_all` ON %s (hfam, org_name, protein_id, go_term_id);" %go_table)
		executeQuery(cursor, "CREATE INDEX `i_org_prot` ON %s (org_name, protein_id);" %go_table)
		db.commit()

	# Retrieving all species from GO table
	logging.info('Checking that all species have GO annotations')
	GOorgs_query = "SELECT DISTINCT org_name FROM %s" %go_table
	(column_names, GOtable_orgs) = executeQuery(cursor, GOorgs_query)
	db.close()

	# Create warning if species are missing InterPro annotation
	GOtable_orgs = map(' '.join, GOtable_orgs)
	missing_orgs = list()
	missing_orgs = list(set(homoTable_orgs)-set(GOtable_orgs))
	if len(missing_orgs) > 0:
		logging.warning("These species do not have GO annotation:")
		for morg in missing_orgs:
			logging.warning("%s" %morg)

	logging.info('Finished %s - runtime :%s' %(go_table, str(datetime.datetime.now()-startTimet_go)))
logging.info('--------------------------------------------------')
logging.info("The program has finished - runtime %s" %str(datetime.datetime.now()-startTimet_1))
logging.info('--------------------------------------------------')

# sys.exit("Exit at line %s\nRuntime %s\n\n\n" %(currentframe().f_lineno, str(datetime.datetime.now()-startTimet_1)))
