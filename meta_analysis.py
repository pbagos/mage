import PythonMeta as PMA
import pandas as pd
import math
import random
import numpy as np
from math import sqrt
import gc
import statistics as st
import itertools
from scipy.stats import norm
from scipy.stats.distributions import chi2
from collections import defaultdict

global genes_p_values, meta_analysis_df, study, num_of_studies, effect_size_list, ind1, ind2, all_genes, means1_table, means2_table, genes_for_ea
genes_p_values = []
meta_analysis_df = []
ind1 = []
ind2 = []
all_genes = []


# With this function we load our data. We need a list which will contain
# the filepaths of each study (txt and tab delimeted)
# and a file which will contain the indices of our control and cases on each study
def split_data(dataframe_list):
    expressions_team1 = []
    expressions_team2 = []

    means1_table = []
    means2_table = []
    for df in dataframe_list:
        # take the unique list annotation symbols
        annot_list = df.loc[1][1:].to_list()
        annot_list_unique = sorted(list(set(annot_list)))

        gene_of_study = df.iloc[2:, 0].reset_index(drop=True)
        gene_of_study = gene_of_study.to_frame()

        # team1
        team_cols1 = list(np.array(np.where(df.loc[1] == annot_list_unique[0]), dtype=int).flatten())
        ind1.append(team_cols1)
        # expressions
        data_team1 = df.iloc[2:][team_cols1].astype(float).reset_index(drop=True)
        # concat expressions with the genes
        new_df1 = pd.concat([gene_of_study, data_team1], axis=1)
        expressions_team1.append(new_df1)

        # team2
        team_cols2 = list(np.array(np.where(df.loc[1] == annot_list_unique[1]), dtype=int).flatten())
        # expressions
        data_team2 = df.iloc[2:][team_cols2].astype(float).reset_index(drop=True)

        ind2.append(team_cols2)
        # concat expressions with the genes
        new_df2 = pd.concat([gene_of_study, data_team2], axis=1)
        expressions_team2.append(new_df2)

        n1 = len(data_team1.columns)
        n2 = len(data_team2.columns)

        y1 = data_team1.mean(axis=1)
        y2 = data_team2.mean(axis=1)

        y1 = pd.DataFrame(y1, columns=['MEANS ΤΕΑΜ 1'])
        y2 = pd.DataFrame(y2, columns=['MEANS ΤΕΑΜ 2'])

        z1 = data_team1.std(axis=1)
        z2 = data_team2.std(axis=1)

        z1 = pd.DataFrame(z1, columns=['STANDARD DEVIATIONS ΤΕΑΜ 1'])
        z2 = pd.DataFrame(z2, columns=['STANDARD DEVIATIONS ΤΕΑΜ 2'])

        n1_arr = [n1] * len(y1)
        n2_arr = [n2] * len(y2)

        n1 = pd.DataFrame(n1_arr, columns=["GROUP SIZE"])
        n2 = pd.DataFrame(n2_arr, columns=["GROUP SIZE"])

        means1_arr = pd.concat([gene_of_study, y1, z1, n1], axis=1)
        means2_arr = pd.concat([y2, z2, n2], axis=1)

        all_genes.append(gene_of_study)

        means1_table.append(means1_arr)
        means2_table.append(means2_arr)
        # calculate the means std_dev and columns
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        # print(es_list)
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)

    return expressions_team1, expressions_team2, means1_table, means2_table


# this function conducts a meta-analysis (Random models, IV-Heg,and SMD)
def calc_metadata(expressions_team2, expressions_team1):
    study = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
             'V', 'W', 'X', 'Y', 'Z']

    all_list2 = []

    for i in range(len(expressions_team1)):
        genes = expressions_team1[i][0]

        x1 = expressions_team1[i].iloc[:, 1:].mean(axis=1)
        x2 = expressions_team1[i].iloc[:, 1:].std(axis=1)

        n1 = len(ind2[i])
        n1 = [n1] * len(x1)

        n1 = pd.DataFrame(n1, columns=["GROUP SIZE"])

        x1 = pd.DataFrame(x1, columns=['MEANS ΤΕΑΜ 1'])
        x2 = pd.DataFrame(x2, columns=['STANDARD DEVIATIONS ΤΕΑΜ 1'])

        all_1 = pd.concat([genes, x1, x2, n1], axis=1)

        y1 = expressions_team2[i].iloc[:, 1:].mean(axis=1)
        y2 = expressions_team2[i].iloc[:, 1:].std(axis=1)

        n2 = len(ind1[i])
        n2 = [n2] * len(y1)

        n2 = pd.DataFrame(n2, columns=["GROUP SIZE"])

        y1 = pd.DataFrame(y1, columns=['MEANS ΤΕΑΜ 2'])
        y2 = pd.DataFrame(y2, columns=['STANDARD DEVIATIONS ΤΕΑΜ 2'])

        all_2 = pd.concat([y1, y2, n2], axis=1)

        all = pd.concat([all_1, all_2], axis=1)
        all_list2.append(all)
    all_info = []
    all_genes2 = []

    for i in all_list2:
        all_info.append(i.values.tolist())
    for i in all_genes:
        all_genes2.append(i.values.tolist())

    items = []

    for j in all_info:
        items.append(j)

    # now everything is in the desired dataFrame,but we need to make the genes as keys in hash tables
    items = list(np.concatenate(items))
    items = np.array(items)

    hash_table = []
    for item in items:
        hash_table.append({item[0]: item[1:]})
    hash_table = np.array(hash_table)

    res = defaultdict(list)
    for sub in hash_table:
        for key in sub:
            res[key].append(sub[key])

    hash_table = dict(res)
    res = {k: v for k, v in hash_table.items() if len(v) >= 2}
    hash_table = res

    x = pd.DataFrame(hash_table.values(), index=hash_table.keys()).T

    text = []
    gene_names = x.columns
    # print(list(gene_names)) # Every column is a gene

    global counter, index

    # print("Gene\tEffect size (Hedge's g)\tStandard_Error\tQ\tI_Squared\tTau_Squared\tp_Q_value\tz_test_value\tp_value\n")
    for i in range(len(x.columns)):
        counter = gene_names[i]
        index = i
        column = x.iloc[:, i].dropna()
        # print(column)
        for i, row in enumerate(column):
            temp = ",".join(str(x) for x in row)
            text.append(study[i] + "," + temp)
            # print(text)

        settings = {"datatype": "CONT",  # for CONTinuous data
                    "models": "Random",  # models: Fixed or Random
                    "algorithm": "IV-Heg",  # algorithm: IV
                    "effect": "SMD"}  # effect size: MD, SMD
        main(text, settings)
        text = []
    return pd.DataFrame(meta_analysis_df)


# This function helps in the execution for the meta_analysis function,and calls the showresults() function
def main(stys, settings):
    d = PMA.Data()  # Load Data class
    m = PMA.Meta()  # Load Meta class
    # f = PMA.Fig()   #Load Fig class

    # You should always tell the datatype first!!!
    d.datatype = settings["datatype"]  # set data type, 'CATE' for binary data or 'CONT' for continuous data
    studies = d.getdata(stys)  # load data
    # studies = d.getdata(d.readfile("studies.txt"))  #get data from a data file, see examples of data files
    # print(showstudies(studies,d.datatype))           #show studies

    m.datatype = d.datatype  # set data type for meta-analysis calculating
    m.models = settings["models"]  # set effect models: 'Fixed' or 'Random'
    m.algorithm = settings["algorithm"]  # set algorithm, based on datatype and effect size
    m.effect = settings["effect"]  # set effect size:RR/OR/RD for binary data; SMD/MD for continuous data
    results = m.meta(studies)  # performing the analysis
    # print(m.models + " " + m.algorithm + " " + m.effect)

    # print (showresults(results))                     #show results table
    showresults(results)
    # f.forest(results).show()                         #show forest plot
    # f.funnel(results).show()                         #show funnel plot


# This function helps in the execution for the meta_analysis function and is called by the main function()
def showresults(rults):
    text = "%-10s %-6s  %-18s %-10s" % ("Study ID", "n", "ES[95% CI]", "Weight(%)\n")

    text1 = str(counter) + "\t" + str(rults[0][1]) + "\t" + str(abs(rults[0][1] / rults[0][10])) + "\t" + str(
        rults[0][7]) + "\t" + str(round(rults[0][9], 2)) + "%\t" + str((rults[0][12])) + "\t" + str(
        (rults[0][8])) + "\t" + str(rults[0][10]) + "\t" + str(rults[0][11])
    new_row = {'Genes': counter, 'p_value': rults[0][11]}
    new_row2 = {'Genes': counter, "Effect size (Hedge's g)": rults[0][1],
                'Standard_Error': abs(rults[0][1] / rults[0][10]), 'Q': rults[0][7], 'I_Squared': round(rults[0][9], 2),
                'Tau_Squared': rults[0][12], 'p_Q_value': rults[0][8], 'z_test_value': rults[0][10],
                'p_value': rults[0][11]}
    genes_p_values.append(new_row)
    meta_analysis_df.append(new_row2)
    # return meta_analysis_df
    return text1


def altmeta(y1, s2):
    n = len(y1)
    w = [(1 / x) for x in s2]
    mu_bar = sum(a * b for a, b in zip(w, y1)) / sum(w)

    Q = sum(a * b for a, b in zip(w, [(x - mu_bar) ** 2 for x in y1]))
    p_Q = chi2.sf(Q, len(y1) - 1)
    H = np.sqrt(Q / (n - 1))
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


# With this function we can conduct a bootstrapped meta_analysis
def bootstrap_analysis(expressions_team1, expressions_team2, means1_table, means2_table, n):
    random_cols = []
    n_of_res = n
    num_of_studies = len(expressions_team1)

    def bootstrap(expressions_team):
        bootstrap_data = []
        for i in range(len(expressions_team)):  # for each study
            g = expressions_team[i][0]  # For each gene
            expr = expressions_team[i].iloc[:, 1:]
            # for each expressions set
            col_list = list(expr.columns)

            for _ in range(n_of_res):  # change it for bigger iterations
                boot_df = pd.DataFrame()

                cols_new = [random.choice(col_list) for _ in range(len(col_list))]
                random_cols.append(cols_new)
                for k in range(len(cols_new)):
                    boot_df = pd.concat([boot_df, expr[cols_new[k]]], axis=1)

                new_df = pd.concat([g, boot_df], axis=1)

                bootstrap_data.append(new_df)

                del boot_df
        return bootstrap_data

    et1 = np.array(bootstrap(expressions_team1), dtype=object).flatten()
    et2 = np.array(bootstrap(expressions_team2), dtype=object).flatten()

    # for step in range(0,len(et1),5):   # change it for bigger iterations
    et1_sliced = [et1[i:i + n_of_res] for i in range(0, len(et1), n_of_res)]
    et2_sliced = [et2[i:i + n_of_res] for i in range(0, len(et2), n_of_res)]

    # print(len(et1_sliced),len(et1_sliced[0]))
    # print((et1_sliced[0][0].std(axis=1)))

    def get_means_and_std(et1_sliced):
        std_boot = []
        means_boot = []
        for i in range(len(et1_sliced)):
            for j in range(len(et1_sliced[0])):
                n = len(et1_sliced[i][j].columns) - 1

                n = [n] * len(et1_sliced[i][j][0])
                n = pd.DataFrame(n, columns=['N'])

                means_boot.append(
                    pd.concat([et1_sliced[i][j][0], et1_sliced[i][j].mean(axis=1), et1_sliced[i][j].std(axis=1),
                               n], axis=1))
        return means_boot, std_boot

    means_boot1, std_boot1 = get_means_and_std(et1_sliced)
    means_boot2, std_boot2 = get_means_and_std(et2_sliced)

    rename_cols1 = ['GENES', 'MEANS1', 'STD1', 'N1']
    rename_cols2 = ['GENES2', 'MEANS2', 'STD2', 'N2']

    boot_data_all = []
    for i in range(len(means_boot1)):
        means_boot1[i].columns = rename_cols1
        means_boot2[i].columns = rename_cols2

    effect_sizes_d = []
    effect_sizes_g = []
    test_g = []
    test_gene = []
    # print(means_boot1)

    for i in range(len(means_boot1)):
        for j in range(len(means_boot1[i])):
            list_of_items1 = list(means_boot1[i].iloc[j])
            list_of_items2 = list(means_boot2[i].iloc[j])
            # print(list_of_items1)
            N = list_of_items1[3] + list_of_items2[3]
            df = N - 2
            J = (math.gamma(df / 2) / (math.sqrt(df / 2) * math.gamma((df - 1) / 2)))
            Sp = sqrt(
                #    (*n1-1)*s1^2                                                (*n2-1)*s2^2
                ((list_of_items1[3] - 1) * (list_of_items1[2] * list_of_items1[2]) + (list_of_items2[3] - 1) * (
                        list_of_items2[2] * list_of_items2[2])) / (list_of_items1[3] + list_of_items2[3] - 2))
            d = (list_of_items1[1] - list_of_items2[1]) / Sp  # effect sizes _d
            g = J * d  # effect size corrected  with Hedge's g

            gc.collect()
            effect_sizes_d.append({list_of_items1[0]: d})
            effect_sizes_g.append({list_of_items1[0]: g})
            test_gene.append(list_of_items1[0])
            test_g.append(g)

    # print(test_g)
    # print(test_gene)
    df = pd.DataFrame(test_g, index=test_gene)
    df = df.sort_index()
    # print(df)
    genes = list(dict.fromkeys(test_gene))  # unique

    temp_df_list = []

    for gene in genes:
        x = df.loc[gene]

        temp_df_list.append(list(x[0]))

    flat_list = list(itertools.chain(*list(n * [np.arange(0, num_of_studies)])))

    # print(flat_list)
    temp_df = pd.DataFrame(temp_df_list, index=genes, columns=flat_list)
    new_cols = np.arange(0, len(list(temp_df.columns)))
    temp_df.columns = new_cols

    # print(temp_df)
    # print(temp_df)
    genes_boot = genes
    # print(genes_boot)
    std_err_boot = []
    for gene in genes_boot:
        tmp = (list(temp_df.loc[gene, :].dropna()))
        tmp_sliced = [tmp[i:i + n_of_res] for i in range(0, len(tmp), n_of_res)]
        # print(tmp_sliced)
        for i in range(len(tmp_sliced)):
            # for j in range (len(tmp_sliced[0])):
            std_err = st.stdev(tmp_sliced[i])
            std_err_boot.append({'Gene': gene, 'std_err': std_err})
    # print(pd.DataFrame(std_err_boot))
    # print(std_err_boot)

    std_err_tmp = pd.DataFrame(std_err_boot, columns=['Gene', 'std_err'])
    std_err_tmp = std_err_tmp.set_index('Gene')

    std_err_list_boot = []
    for gene in genes:
        x = std_err_tmp.loc[gene]

        std_err_list_boot.append(list(x.values.flatten()))

    flat_list = list(itertools.chain(*list([np.arange(0, num_of_studies)])))

    df = pd.DataFrame(std_err_list_boot, index=genes, columns=flat_list)

    boot_df = pd.DataFrame()
    for gene in genes:
        # print(gene)
        # print(df.T.loc[gene,:].dropna())
        tmp = list(df.loc[gene, :].dropna())
        new_row = {'GENES': gene, "SE": tmp}
        # append row to the dataframe
        boot_df = boot_df.append(new_row, ignore_index=True)  # Boot DF has the standard_errors

    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)

    genes_uniq = (boot_df['GENES'])

    es_df = pd.DataFrame()

    for i in range(len(means1_table)):
        for gene in means1_table[i].index:
            # print(list(means2_table[i].iloc[gene]))
            m1 = (list(means1_table[i].iloc[gene])[1])
            m2 = (list(means2_table[i].iloc[gene])[0])

            # print(list_of_items1)
            n1 = list(means1_table[i].iloc[gene])[3]
            n2 = list(means2_table[i].iloc[gene])[2]
            N = n1 + n2
            df = N - 2
            J = (math.gamma(df / 2) / (math.sqrt(df / 2) * math.gamma((df - 1) / 2)))
            Sp = sqrt(
                #    (*n1-1)*s1^2                                                (*n2-1)*s2^2
                ((n1 - 1) * (list(means1_table[i].iloc[gene])[2] * list(means1_table[i].iloc[gene])[2]) + (n2 - 1) * (
                        list(means2_table[i].iloc[gene])[1] * list(means2_table[i].iloc[gene])[1])) / (N - 2))
            d = (m1 - m2) / Sp  # effect sizes _d
            g = J * d  # effect size corrected  with Hedge's g
            var_d = N / (n1 * n2) + (d * d) / (2 * N)
            var_g = J * J * var_d
            new_row = {'GENES': (list(means1_table[i].iloc[gene])[0]), 'd': d, 'g': g, 'var_d': var_d, 'var_g': var_g}
            es_df = es_df.append(new_row, ignore_index=True)
    # print(es_df)
    es_list = pd.DataFrame()
    for gene in genes_uniq:
        es_df_index = es_df[es_df['GENES'] == gene]
        es_of_studies = list(es_df_index['g'])  # or d
        new_row = {"GENES": gene, "List_of_ES": es_of_studies}
        es_list = es_list.append(new_row, ignore_index=True)
        # print(es_df_index)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    # print(es_list)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    # print(boot_df)

    total_df = pd.DataFrame()

    list_of_boot = []
    # for each unique gene
    # Take each line  from the  to es_list and from  boot_df and do the meta-analysis
    # print('Genes'+'\t'+"Effect Size (Hedge's g)"+'\t'+'Standard Error '+'\t'+'Q'+'\t'+'I - Squared'+'\t'+'Tau - Squared'+'\t'+'p_Q'+'\t'+'\t'+'z-score'+'\t'+'p-value')
    num_of_chunks = 10
    es_list_chunks = [es_list[i:i + num_of_chunks] for i in range(0, len(es_list), num_of_chunks)]
    boot_df_chunks = [boot_df[i:i + num_of_chunks] for i in range(0, len(boot_df), num_of_chunks)]
    for es_list, bood_df in zip(es_list_chunks, boot_df_chunks):
        for index in es_list.index:
            gene = (es_list.loc[index, 'GENES'])
            d = (es_list.loc[index, 'List_of_ES'])
            se = boot_df.loc[index, 'SE']
            # print(d)
            #
            # print(se)

            if (len(se) > 1):
                # Gene	Effect size (Hedge's g)	Q	I - Squared	Tau - Squared	p_Q value	z-test value	p-value
                # Genes'+'\t'+"Effect Size"+'\t'+'Q'+'\t'+'p_Q'+'\t'+'I - Squared'+'\t'+'z-score'+'\t'+'p-value'

                Q, I2, tau2, p_Q, st_e, z, e_size, p = altmeta(d, se)

                # print(gene+"\t"+str(e_size)+"\t"+str(st_e)+"\t"+str(Q)+"\t"+str(I2)+"\t"+str(tau2)+"\t"+str(p_Q)+"\t"+str(z)+"\t"+str(p))
                #
                new_row2 = {'Genes': gene, "Effect size (Hedge's g)": e_size, 'Standard_Error': st_e, 'Q': Q,
                            'I_Squared': I2, 'Tau_Squared': tau2, 'p_Q_value': p_Q, 'z_test_value': z, 'p_value': p}
                list_of_boot.append(new_row2)

    # print(total_df)
    return pd.DataFrame(list_of_boot)


# With this function we can get the  one step methods (Bonferroni and Sidak)
def get_one_step_methods(meta_analysis_df,
                         alpha):  # This function takes the Sidak and the Bonferoni multiple test method

    m = len(meta_analysis_df['p_value'])
    p_values = (np.array(meta_analysis_df['p_value'], dtype=float))
    bonf = []
    sidak = []
    exponent = 1 / m
    for i in range(m):  # a = 0.05 due to the 95 % CI
        if p_values[i] >= (1 - (1 - alpha) ** exponent):
            sidak.append(0)  # 0 ---> no significant difference
        else:
            sidak.append(1)  # 1 --->   significant  difference

        if p_values[i] >= alpha / m:
            bonf.append(0)

        else:
            bonf.append(1)

    one_step_methods = pd.DataFrame(list(zip(list(meta_analysis_df['Genes']), p_values, bonf, sidak)),
                                    columns=['genes_one_step', 'p_values_one_step', 'bonferroni', 'sidak'])
    # print(one_step_methods)
    return (one_step_methods)


# With this function we can get the  step down methods (Holm and Holland)
def get_step_down_methods(meta_analysis_df,
                          alpha):  # This function takes the Sidak and the Bonferoni multiple test method

    m = len(meta_analysis_df['p_value'])
    meta_analysis_df = meta_analysis_df.sort_values(by=['p_value'])
    # print(genes_and_p_values)
    p_values = (np.array(meta_analysis_df['p_value'], dtype=float))

    holm = []
    holland = []
    for i in range(m):
        if p_values[i] >= (alpha / (m - (i) + 1)):

            holm.append(0)
        else:
            holm.append(1)

        exponent = 1 / (m - (i) + 1)
        if p_values[i] >= (1 - (1 - alpha) ** exponent):

            holland.append(0)

        else:
            holland.append(1)

    step_down_methods = pd.DataFrame(list(zip(list(meta_analysis_df['Genes']), p_values, holm, holland)),
                                     columns=['genes_step_down', 'p_values_step_down', 'holm', 'holland'])
    # print(step_down_methods)
    return (step_down_methods)


# With this function we can get the step up methods (Simes and Hochberg)
def get_step_up_methods(meta_analysis_df,
                        alpha):  # This function takes the Sidak and the Bonferoni multiple test method

    m = len(meta_analysis_df['p_value'])
    meta_analysis_df = meta_analysis_df.sort_values(by=['p_value'], ascending=False)  # sort by p_value
    # print(genes_and_p_values)
    p_values = (np.array(meta_analysis_df['p_value'], dtype=float))

    hochberg = []
    simes = []
    for i in range(m):

        if p_values[i] >= (alpha / (m - (i) + 1)):
            hochberg.append(0)
        else:
            hochberg.append(1)

        if p_values[i] >= ((i * alpha) / m):
            simes.append(0)
        else:
            simes.append(1)

    step_down_methods = pd.DataFrame(list(zip(list(meta_analysis_df['Genes']), p_values, hochberg, simes)),
                                     columns=['genes_step_down', 'p_values_step_down', 'hochberg', 'simes'])
    # print(step_down_methods)
    return step_down_methods


# call  all the Multiple - tests functions (Bonferroni,Sidak,Holm,Holland,Simes and Hochberg from the upper functions)
def all_multiple_tests(meta_analysis_df, alpha):
    d1 = get_one_step_methods(meta_analysis_df, alpha)
    d2 = get_step_down_methods(meta_analysis_df, alpha)
    d3 = get_step_up_methods(meta_analysis_df, alpha)

    total_df = pd.concat([d1, d2, d3], axis=1)
    return total_df


def run(settings, data):
    num_of_reps = int(settings['num_of_reps'])
    alpha = float(settings['alpha'])
    mult_tests = settings['multiple_comparisons']
    bootstrap = settings['bootstrap']

    # Splits the cases and the controls of our study
    expressions_team1, expressions_team2, means1_table, means2_table = split_data(data)

    if bootstrap == 'YES':
        print("Bootstrap Option")
        meta_analysis_df = bootstrap_analysis(expressions_team2, expressions_team1, means1_table, means2_table,
                                              n=num_of_reps)
    else:
        meta_analysis_df = calc_metadata(expressions_team1, expressions_team2)

    # Each function of these functions, conducts the multiple test functions  and returns a dataframe.
    step_down = get_step_down_methods(meta_analysis_df, alpha)
    step_up = get_step_up_methods(meta_analysis_df, alpha)
    one_step = get_one_step_methods(meta_analysis_df, alpha)

    # #   all multiple tests  (bonferroni, sidak, holm, holland, hochberg and simes methods)
    # #   0 --> No statistical significance difference
    # #   1 --> Statistical significance difference

    all_tests = all_multiple_tests(meta_analysis_df, alpha)

    if mult_tests == "one_step":
        tests = one_step
    elif mult_tests == "step_up":
        tests = step_up
    elif mult_tests == "step_down":
        tests = step_down
    else:
        tests = all_tests

    meta_an = pd.concat([meta_analysis_df, tests], axis=1)

    return meta_an
