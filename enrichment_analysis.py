import requests
import pandas as pd
import numpy as np
from gprofiler import GProfiler


def run(settings, genes):
    p_organism = settings['organism']
    p_threshold = float(settings['threshold'])
    p_threshold_method = settings['threshold_method']

    # genes = ['VKORC1', 'ABCB1', 'ABCC2', 'C18orf56', 'CLCN6', 'CYP2A7P1', 'CYP2B6', 'CYP2C8', 'CYP3A4', 'CYP3A5',
    #          'DRD3', 'FCGR2A', 'FCGR3A', 'G6PD', 'HLA-B', 'HLA-DRB1', 'KCNJ2', 'MTHFR', 'MTRR', 'NME4', 'NR1I2',
    #          'SLCO1B1', 'SOD1', 'SORCS2', 'TYMS', 'UGT1A1', 'UGT1A10', 'UGT1A3', 'UGT1A4', 'UGT1A5', 'UGT1A6', 'UGT1A7',
    #          'UGT1A8', 'UGT1A9', 'ZSCAN25']

    #genes =['ABCB1','ABCC1','ABCC2','APOC1','APOC3','APOE','B4GALT2','C18orf56','TYMSOS','CTLA4','CYP2C19','CYP2C9','CYP2D6','CYP3A4','CYP3A5','CYP4F2','ENOSF1','G6PD','GPX1','HLA-DQA1','HLA-DRB1','IRS1','ITGA2','ITGB3','MT-ND3','MTR','NAT2','NOS3','NR1I2','NTRK1','PEAR1','PTGS1','SLC19A1','SLCO1B1','SLCO2B1','SORCS2','TOMM40','TYMS','UGT1A1','UGT1A3','UGT1A6','UGT1A7','VDR','VEGFA','ZSCAN25']

    gp = GProfiler(return_dataframe=True)
    var = gp.profile(organism=p_organism,
                     query=genes,
                     all_results=True,
                     ordered=True,
                     combined=False,
                     measure_underrepresentation=False,
                     no_iea=False,
                     domain_scope='known',
                     numeric_namespace='AFFY_HUGENE_1_0_ST_V1',
                     sources=["GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "WP", "TF", "MIRNA", "HPA", "CORUM", "HP"],
                     user_threshold=p_threshold,
                     significance_threshold_method=p_threshold_method,
                     no_evidences=False
                     )

    var['negative_log10_of_adjusted_p_value'] = -np.log10(var.p_value)
    #print(var.head().to_string())
    #print(var.head())
    return var
