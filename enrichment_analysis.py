import numpy as np
from gprofiler import GProfiler


def run(settings, genes):
    p_organism = settings['organism']
    p_threshold = float(settings['threshold'])
    p_threshold_method = settings['threshold_method']

    gp = GProfiler(return_dataframe=True)
    var = gp.profile(organism=p_organism,
                     query=genes,
                     all_results=True,
                     ordered=False,
                     combined=False,
                     measure_underrepresentation=False,
                     no_iea=False,
                     no_evidences=False,
                     domain_scope='annotated',
                     numeric_namespace='AFFY_HUGENE_1_0_ST_V1',
                     sources=["GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "WP", "TF", "MIRNA", "HPA", "CORUM", "HP"],
                     user_threshold=p_threshold,
                     significance_threshold_method=p_threshold_method,
                     )

    var['negative_log10_of_adjusted_p_value'] = -np.log10(var.p_value)

    #print(var.head())
    return var
