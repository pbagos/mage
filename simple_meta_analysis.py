import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.stats.distributions import chi2
from meta_analysis import all_multiple_tests, get_one_step_methods, get_step_down_methods, get_step_up_methods


def altmeta(y1, s2):
    n = len(y1)
    w = [(1 / x) for x in s2]
    mu_bar = sum(a * b for a, b in zip(w, y1)) / sum(w)

    Q = sum(a * b for a, b in zip(w, [(x - mu_bar) ** 2 for x in y1]))
    p_Q = chi2.sf(Q, len(y1) - 1)
    H = np.sqrt(Q / (n - 1))
    if (Q == 0):
        I2 = 0
        tau2_DL = 0
    else:
        I2 = (Q - (len(s2) - 1)) / Q
        I2 = max(0, I2)
        tau2_DL = (Q - n + 1) / (sum(w) - sum([x ** 2 for x in w]) / sum(w))
        tau2_DL = max(0, tau2_DL)

    # re adjustment of the weights
    w = [(1 / (x + tau2_DL)) for x in s2]
    mu_bar = sum(a * b for a, b in zip(w, y1)) / sum(w)
    se = np.sqrt(1 / sum(w))

    z = (mu_bar / se)
    p = 1 - norm.cdf(abs(z))
    return (Q, I2, tau2_DL, p_Q, se, z, mu_bar, p)


def simple_meta_analysis(file_list, folder_path, alpha, mult_tests):
    # Read each file as a pandas DataFrame with tab delimiter
    dataframes = []
    results_cols = ['Genes', "Effect size (Hedge's g)", 'Standard_Error', 'Q', 'I_Squared', 'Tau_Squared', 'p_Q_value',
                    'z_test_value', 'p_value', 'num_of_studies']
    results_list = []

    for file_path in file_list:
        df = pd.read_csv(folder_path + file_path, delimiter='\t', header=None)

        # Add the DataFrame to the list
        dataframes.append(df)

    # Combine all DataFrames into a single DataFrame
    combined_df = pd.concat(dataframes, ignore_index=True)

    # From the combined dataframe fetch: Gene, Effect_size, Standard error and n1,n2
    # and perform altmeta for each gene

    # Convert dataframe to dictionary

    gene_dict = {}
    for row in combined_df.itertuples(index=False):
        gene_id = row[0]
        values = list(row[1:])
        if gene_id in gene_dict:
            gene_dict[gene_id].append(values)
        else:
            gene_dict[gene_id] = [values]
    del gene_dict['Gene']

    # print(gene_dict)

    for gene in gene_dict.keys():
        data = pd.DataFrame(gene_dict[gene])
        es_list = list(data[0])
        se_list = list(data[1])

        # Convert to floats
        es_list = [float(num) for num in es_list]
        se_list = [float(num) for num in se_list]

        Q, I2, tau2_DL, p_Q, se, z, mu_bar, p = altmeta(es_list, se_list)
        # ['Genes', "Effect size (Hedge's g)", 'Standard_Error', 'Q', 'I_Squared', 'Tau_Squared', 'p_Q_value',
        #  'z_test_value', 'p_value', 'num_of_studies', 'cases', 'controls']
        results_list.append([gene, mu_bar, se, Q, I2, tau2_DL, p_Q, z, p, len(es_list)])

    results_df = pd.DataFrame(results_list, columns=results_cols)

    step_down = get_step_down_methods(results_df, alpha)
    step_up = get_step_up_methods(results_df, alpha)
    one_step = get_one_step_methods(results_df, alpha)

    all_tests = all_multiple_tests(results_df, alpha)

    if mult_tests == "one_step":
        tests = one_step
    elif mult_tests == "step_up":
        tests = step_up
    elif mult_tests == "step_down":
        tests = step_down
    else:
        tests = all_tests

    meta_an = pd.concat([results_df, tests], axis=1)

    meta_an = meta_an.drop(
        ['genes_one_step', 'genes_step_down', 'genes_step_up', 'p_values_step_down', 'p_values_step_up',
         'p_values_one_step'], axis=1)

    return meta_an
