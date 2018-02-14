"""
Download info on secondary metabolic gene clusters using jgi codes.
"""


import os
import csv
import sys
import pandas as pd
from bioSlim3 import dbFetch, dbwHeader




def dlSMdata(orgSet):
    """
    Download smurf gff and interpro data and combine them.
    orgSet has to be a list of jgi names.
    """
    print("Combining smurf, gff and organism data")
    query = """
    SELECT org.org_id, org.name, org.genus, org.real_name, org.section,
    smurf.sm_protein_id AS protein_id, smurf.sm_short, CONCAT(smurf.org_id, '_', smurf.clust_backbone, '_', smurf.clust_size) AS cluster_id, smurf.clust_size,
    gff_prot.gff_start, gff_prot.gff_end, gff_prot.gff_strand
    FROM smurf
    LEFT JOIN (SELECT MIN(gff.gff_start) AS gff_start, MAX(gff.gff_end) AS gff_end, gff.org_id, gff.gff_protein_id, gff.gff_strand, gff.gff_seqorigin FROM gff GROUP BY gff.org_id, gff.gff_protein_id) AS gff_prot ON
    smurf.org_id = gff_prot.org_id AND smurf.sm_protein_id = gff_prot.gff_protein_id
    JOIN organism AS org ON smurf.org_id = org.org_id AND org.name IN ('%s');
    """ % "','".join(orgSet)

    smurfGff = dbwHeader(query)

    print("Downloading interpro annotations")

    query = """
    SELECT phi.org_id, phi.protein_id,
    GROUP_CONCAT(DISTINCT ipr.ipr_desc ORDER BY phi.ipr_domain_start SEPARATOR ',') AS ipr_desc ,
    GROUP_CONCAT(DISTINCT ipr.ipr_domaindesc ORDER BY phi.ipr_domain_start SEPARATOR ',') AS ipr_domaindesc , ipr.ipr_domaindb
    FROM protein_has_ipr AS phi
    JOIN ipr ON phi.ipr_id = ipr.ipr_id AND phi.ipr_domaindb = 'HMMPfam' AND phi.ipr_domaindb = ipr.ipr_domaindb
    GROUP BY phi.org_id, phi.protein_id;
    """

    ipr = dbwHeader(query)

    smurf = pd.DataFrame(smurfGff[1], columns = smurfGff[0])
    ipr = pd.DataFrame(ipr[1], columns = ipr[0])

    print("Merging data")
    smIpr = pd.merge(smurf, ipr,  how='left', left_on=['org_id','protein_id'], right_on = ['org_id','protein_id'])

    return(smIpr)
