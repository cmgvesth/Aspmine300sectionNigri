#!/usr/bin/python3

"""
Main secondary metabolism pipeline. Imports modules and processes data to write a dataframe which can later be processed with R.
"""

import os
import csv
import sys
import argparse
import pandas as pd
import shutil
sys.path.append(os.path.join(os.getcwd(),"scripts"))

from smServerSide import tmpSmBiTable, createBidirSmurf, mysqlSmChecker
from aspSMDl import dlSMdata
from bioSlim3 import tupleToFasta
from processMibig3 import processMibig, dlSmurfProteins, writeMibigFormatted, processBlastResult, dlMibig

import misc



parser=argparse.ArgumentParser(description='''
Script to execute the seconary metabolite analysis pipeline. In case you want to leave out some analysis you have to modify the script.\n

Example:

python3 smDataHandling.py -o flavi_orgs.txt -bibase I_Flavi_biblast -biFinal smurf_bidir_hits_flavi -t ingek_tree_JGIname.nwk -l flavi.log -s flavi_test
''')
parser.add_argument("--orgs", "-o",
					dest="filename",
					required=True,
					help="Input file with jgi names (one per row) of organisms", metavar="FILE")
parser.add_argument("--treeFile", "-t",
					dest="tree",
					required=True,
					help="Tree file in newick format",
					metavar="FILE")
parser.add_argument("--biblastTable", "-bibase",
					dest="bibase",
					required=True,
					help="Specify the original biblast table to use as basis for the smurf bidirectional hits table", metavar="CHAR")
parser.add_argument("--clusterBiblast", "-biFinal",
					dest="biFinal",
					required=True,
					help="Specify name for smurf bidirectional hits table", metavar="CHAR")
parser.add_argument("--log", "-l",
					dest="logFile",
					required=True,
					help="Specify the name for log file", metavar="CHAR")
parser.add_argument("--outdir", "-od",
					dest="sn",
					required=True,
					help="Output directory, preferably the name of your set", metavar="CHAR")
parser.add_argument("--legacyBlast", "-lb",
					dest="legacyBlast",
					action="store_true",
					default=False,
					help="Specify the original biblast table to use as basis for the smurf bidirectional hits table")
args=parser.parse_args()




with open("config.txt") as c:
	config = misc.readConfig(c.readlines())

print("Queries to format smurf databases")
if args.legacyBlast:
	print("%s -i smurf.fasta" % config['formatDbPath'])
else:
	print("%s -in smurf.fasta -dbtype prot" % config['makeblastdbPath'])


########
# LOADING ORGS

filename = args.filename

biblastBaseTable = args.bibase
smurfBidirHitsName = args.biFinal
testLogName = args.logFile
treeFile = args.tree
setName = args.sn

with open(filename, "r") as tmp:
	orgSet = [item.strip() for item in tmp.readlines()]

if setName not in os.listdir():
	os.mkdir(setName)

else:
	input("A folder for the specified set is already available. If you want to rerun the analysis press Enter, else ctrl+c/d\n")

shutil.copy2(treeFile, setName)

os.chdir(setName)
# Checking data
mysqlSmChecker(orgSet, testLogName)

# Creating smurf bidir hits table for dataset.
tmpSmBiTable(smtable = biblastBaseTable) # Creating a temporary table

createBidirSmurf(smtable = smurfBidirHitsName) # Creating a sm cluster table


# DOWNLOADING SM DATA
input("Did you create a smurfbiblast table? If yes, press Enter, if not, ctrl+c/d")
smIpGf = dlSMdata(orgSet)






##########
# MIBIG PART

input("Press Enter to continue with dereplication of mibig entries in your dataset")

print("Creating directory for mibig")

if "mibig" not in os.listdir():
	os.mkdir("mibig")

os.chdir("mibig")
dlMibig()

sp = dlSmurfProteins(orgSet)
tupleToFasta(sp,"smurf.fasta")
print("Wrote smurf proteins to disk")

try:
	print("Formatting smurf database")
	if args.legacyBlast:
		os.system("%s -i smurf.fasta" % config['formatDbPath'])
	else:
		os.system("%s -in smurf.fasta -dbtype prot" % config['makeblastdbPath'])
except Exception as e:
	print("Could not format smurf proteins to blast database")

print(os.getcwd())
print("processing mibig data")
geneClusters = processMibig("mibigGbkFiles")
print("Writing mibig data to disk")
writeMibigFormatted(geneClusters)


input("""Press Enter to perform blast: \n
%s -query mibigDb.fasta -db smurf.fasta -max_target_seqs 25 -evalue 1e-50 -num_threads 3 -outfmt '6 std qlen slen' -out mibigVsSmurf.txt\n WILL TAKE APPROX 30MIN TO 1H.""" % config['blastpPath'])
# os.getcwd()
os.system("%s -query mibigDb.fasta -db smurf.fasta -max_target_seqs 25 -evalue 1e-50 -num_threads 3 -outfmt '6 std qlen slen' -out mibigVsSmurf.txt" % config['blastpPath'])

print("Processing mibig blast results")
mibigBlast = processBlastResult(blastFile = "mibigVsSmurf.txt", translationFile = "translateBgc.txt")

mibigBlast['org_id'] = mibigBlast['org_id'].astype('int64')

mibigBlast['protein_id'] = mibigBlast['protein_id'].astype('int64')

os.chdir("..")


try:
    print("Joining mibig annotation on sm data")
    smData = pd.merge(smIpGf, mibigBlast[['org_id', 'protein_id', 'compound']], how = "left", left_on=['org_id','protein_id'], right_on = ['org_id','protein_id'])
except Exception as e:
    raise
    print("Cannot merge mibig results and sm data")

smData.to_csv("sm_data_"+setName+".tsv", sep = '\t', encoding = "UTF-8", index = False, na_rep='none')


# Using r scipt to create cluster families
os.chdir("..")

print("Current working directory")
print(os.getcwd())

smFile = "sm_data_"+setName+".tsv"

cmdString = "Rscript clusterData.R %s %s %s" % (smFile, setName, smurfBidirHitsName) # File was clustered before, so it will be renamed with a _c extension
print(cmdString)
print("Calculating cluster families")

os.system(cmdString)



print("Searching for unique secondary metabolic gene clusters in tree")
cmdString = "Rscript uniquesAtNodes.R %s %s %s" % (smFile.replace(".tsv","")+"_c.tsv", setName, treeFile)

os.system(cmdString)
