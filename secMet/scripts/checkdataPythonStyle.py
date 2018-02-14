#!/usr/bin/python3

import sys
import logging
import bioSlim3 as bio
import argparse




def parseOrgProt(l):
    d = {}
    for org, prot in l:
        if str(org) not in d:
            d[str(org)] = []
        d[str(org)].append(str(prot))
    return(d)

# https://fangpenlin.com/posts/2012/08/26/good-logging-practice-in-python/

def mysqlSmChecker(orgSet, logfile):
    logging.basicConfig(filename=logfile,level=logging.DEBUG)

    # Getting org_ids from organism table
    handle = bio.dbFetch("""SELECT org_id, name, real_name, section  FROM organism WHERE name IN ('%s')""" % "','".join(orgSet))

    orgIds = {}
    for org_id, name, real_name, section in handle:
        try:
            int(org_id)
            orgIds[name] = str(org_id)
        except Exception as e:
            print("Species %s does not have an valid org_id" % name)
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

    # IPR checking

    print("Downloading ipr data")
    ipr_idsPerOrgRaw = bio.dbFetch("""SELECT org_id, ipr_id FROM protein_has_ipr WHERE org_id IN ('%s') GROUP BY org_id, ipr_id""" % "','".join(org_ids) )

    ipr_idsPerOrg = parseOrgProt(protein_has_iprRaw)

    iprRaw = bio.dbFetch("""SELECT DISTINCT(ipr_id) FROM ipr """)
    ipr = [item[0] for item in iprRaw]

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
                    logging.warning("%s Duplicate sm_protein_ids for the following cluster backbones: %s" % (name,",".join(multiBbPerProt)))
                    print("%s Duplicate sm_protein_ids for the following cluster backbones: %s" % (name,",".join(multiBbPerProt)))

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
            try:
                phiCheck = set(ipr) - set(ipr_idsPerOrg[org_id])
                if gffCheck:
                    logging.warning("Missing ipr entries for %s", name)
            except Exception as e:
                logging.error("Could not get ipr files for %s", name)
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

# if __name__ == '__main__':
#     parser=argparse.ArgumentParser(
#         description='''Script checks smurf, proteins and protein_has_ipr tables for shared protein ids.
#         Results are written to checkSmDb.log.
#         Uses python 3.''')
#     args=parser.parse_args()
#     with open("../testSet.txt") as testf:
#         orgSet = [line.strip() for line in testf]
