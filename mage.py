import argparse
import configparser
import pandas as pd
import numpy as np
import time

# Custom Libraries
import gisu
import meta_analysis
import multivariate
import plots
import enrichment_analysis

global settings, gprofiler_settings, version, studies
settings = {}
studies = []
studies_transform = []
version = '1.0.2'


def parse_args():
    print("Preparing System Arguments");
    """Parse input arguments."""
    parser = argparse.ArgumentParser(description='Process Tool Arguments')
    # One argument for the configuration file conf.txt
    parser.add_argument('-c', metavar='--conf', required=True, help='Configuration File', type=str, default='conf.txt')
    parser.add_argument('-o', metavar='--output', required=True, help='Output Directory', type=str, default='results/')

    args = parser.parse_args()
    return args


def parse_conf(conf_filename):
    print("Preparing System Configuration (" + conf_filename + ")");
    """Parse configuration arguments."""
    conf = configparser.ConfigParser(comment_prefixes='#', allow_no_value=True)
    conf.read(conf_filename)

    for key, val in conf.items('SETTINGS'):
        settings[key] = val


if __name__ == '__main__':
    t0 = time.time()
    print("MAGE :: Meta-Analysis of Gene Expression")
    print("Version " + version + "; September 2021")
    print("Copyright (C) 2021 Pantelis Bagos")
    print("Freely distributed under the GNU General Public Licence (GPLv3)")
    print("--------------------------------------------------------------------------")
    # initialization step
    args = parse_args()
    parse_conf(args.c)
    filepath = args.o

    print("Load File data")
    file_list = list(settings['study_files'].split(","))

    for i in range(len(file_list)):
        # Read file data
        studypath = settings['study_dir'] + file_list[i].strip()
        file = pd.read_csv(studypath, sep='\t', low_memory=False, header=None)
        studies.append(file)

    if settings.get('run_gisu') == 'YES':
        print("Gene ID/Symbol update started")
        platforms = list(settings['platform'].split(","))

        for study in studies:
            study_transform = gisu.run(settings, study, platforms[i])
            studies_transform.append(study_transform)
        data = studies_transform
    else:
        data = studies

    if settings['multivariate'] == 'YES':
        print('Multivariate Analysis started')
        metanalysis_df = multivariate.run(settings, data, filepath)
        metanalysis_df.to_csv(filepath + 'multivariate_analysis_results.txt', sep='\t', mode='w')
    else:
        print('Meta Analysis started')
        metanalysis_df = meta_analysis.run(settings, data)
        metanalysis_df.to_csv(filepath + 'meta_analysis_results.txt', sep='\t', mode='w')

        # create and save plots
        if settings.get('plots') == 'YES':
            plots.meta_analysis_plots(metanalysis_df, filepath)

    if settings.get('enrichment_analysis') == 'YES':
        print('Enrichment Analysis started')
        alpha = float(settings['alpha'])
        genes_for_ea = metanalysis_df['Genes'].where(
            np.array(metanalysis_df['p_value'], dtype=float) < alpha).dropna().tolist()

        enrichment_analysis_df = enrichment_analysis.run(settings, genes_for_ea)
        enrichment_analysis_df.to_csv(filepath + 'enrichment_analysis_results.txt',
                                      header=enrichment_analysis_df.columns, index=None, sep='\t', mode='a')
        if settings.get('plots') == 'YES':
            plots.ea_manhattan_plot(enrichment_analysis_df, filepath)
            plots.ea_heatmap_plot(enrichment_analysis_df, filepath)

    t1 = time.time()
    total = t1 - t0
    print("Execution time = " + str(total) + '\t' + '  seconds')