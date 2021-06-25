import pandas as pd
import math
import scipy.stats as st
import numpy as np
import plots

from scipy.stats import norm
from scipy.stats import chi2
from collections import defaultdict

global ind1, ind2, ind3, all_genes, num_of_studies

all_genes = []
ind1 = []
ind2 = []
ind3 = []


# With this function we can get the  one step methods (Bonferroni and Sidak)
def get_one_step_methods(meta_analysis_df,
                         alpha, name):  # This function takes the Sidak and the Bonferoni multiple test method
    m = len(meta_analysis_df[name])
    p_values = (np.array(meta_analysis_df[name], dtype=float))
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
                                    columns=['genes_one_step', name, 'bonferroni_' + name, 'sidak_' + name])

    return (one_step_methods)


# With this function we can get the  step down methods (Holm and Holland)
def get_step_down_methods(meta_analysis_df,
                          alpha, name):  # This function takes the Sidak and the Bonferoni multiple test method
    m = len(meta_analysis_df[name])
    meta_analysis_df = meta_analysis_df.sort_values(by=[name])
    p_values = (np.array(meta_analysis_df[name], dtype=float))

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
                                     columns=['genes_step_down', name, 'holm_' + name, 'holland_' + name])

    return (step_down_methods)


# With this function we can get the step up methods (Simes and Hochberg)
def get_step_up_methods(meta_analysis_df,
                        alpha, name):  # This function takes the Sidak and the Bonferoni multiple test method
    m = len(meta_analysis_df['p_value'])
    meta_analysis_df = meta_analysis_df.sort_values(by=[name], ascending=False)  # sort by p_value
    # print(genes_and_p_values)
    p_values = (np.array(meta_analysis_df[name], dtype=float))

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
                                     columns=['genes_step_up', name, 'hochberg_' + name, 'simes_' + name])
    # print(step_down_methods)
    return (step_down_methods)


def altmeta(y1, s2):
    n = len(y1)
    w = [(1 / x) for x in s2]

    mu_bar = sum(a * b for a, b in zip(w, y1)) / sum(w)

    Q = sum(a * b for a, b in zip(w, [(x - mu_bar) ** 2 for x in y1]))

    tau2_DL = (Q - n + 1) / (sum(w) - sum([x ** 2 for x in w]) / sum(w))
    tau2_DL = max(0, tau2_DL)

    # re adjustment of the weights
    w = [(1 / (x + tau2_DL)) for x in s2]
    mu_bar = sum(a * b for a, b in zip(w, y1)) / sum(w)
    se = np.sqrt(1 / sum(w))

    z = (mu_bar / se)
    p = 2 * (1 - norm.cdf(abs(z)))

    return (se, mu_bar, p)


def split_data(dataframe_list):
    expressions_team1 = []
    expressions_team2 = []
    expressions_team3 = []

    means1_table = []
    means2_table = []
    means3_table = []

    all_genes = []
    for df in dataframe_list:
        # take the unique list annotation symbols
        annot_list = df.loc[1][1:].to_list()
        annot_list_unique = sorted(list(set(annot_list)))

        gene_of_study = df.iloc[2:][0].reset_index(drop=True)
        gene_of_study = gene_of_study.to_frame()

        gene_of_study2 = gene_of_study.rename(columns={0: 'GENES'}, inplace=False)

        # team1
        team_cols1 = list(np.array(np.where(df.loc[1] == annot_list_unique[0]), dtype=int).flatten())
        ind1.append(team_cols1)
        # expressions
        data_team1 = df.iloc[2:][team_cols1].astype(float).reset_index(drop=True)
        # concat expressions with the genes
        new_df1 = pd.concat([gene_of_study, data_team1], axis=1)
        expressions_team1.append(new_df1)

        # calculate the means std_dev and columns

        # team2
        team_cols2 = list(np.array(np.where(df.loc[1] == annot_list_unique[1]), dtype=int).flatten())
        # expressions
        data_team2 = df.iloc[2:][team_cols2].astype(float).reset_index(drop=True)

        ind2.append(team_cols2)
        # concat expressions with the genes
        new_df2 = pd.concat([gene_of_study, data_team2], axis=1)
        expressions_team2.append(new_df2)

        # team3
        team_cols3 = list(np.array(np.where(df.loc[1] == annot_list_unique[2]), dtype=int).flatten())
        # expressions
        data_team3 = df.iloc[2:][team_cols3].astype(float).reset_index(drop=True)

        ind3.append(team_cols3)
        # concat expressions with the genes
        new_df3 = pd.concat([gene_of_study, data_team3], axis=1)
        expressions_team3.append(new_df3)

        n1 = len(data_team1.columns)
        n2 = len(data_team2.columns)
        n3 = len(data_team3.columns)

        y1 = data_team1.mean(axis=1)
        y2 = data_team2.mean(axis=1)
        y3 = data_team3.mean(axis=1)

        y1 = pd.DataFrame(y1, columns=['MEANS ΤΕΑΜ 1'])
        y2 = pd.DataFrame(y2, columns=['MEANS ΤΕΑΜ 2'])
        y3 = pd.DataFrame(y3, columns=['MEANS ΤΕΑΜ 3'])

        z1 = data_team1.std(axis=1)
        z2 = data_team2.std(axis=1)
        z3 = data_team3.std(axis=1)

        z1 = pd.DataFrame(z1, columns=['STANDARD DEVIATIONS ΤΕΑΜ 1'])
        z2 = pd.DataFrame(z2, columns=['STANDARD DEVIATIONS ΤΕΑΜ 2'])
        z3 = pd.DataFrame(z3, columns=['STANDARD DEVIATIONS ΤΕΑΜ 3'])

        n1_arr = [n1] * len(y1)
        n2_arr = [n2] * len(y2)
        n3_arr = [n3] * len(y2)

        n1 = pd.DataFrame(n1_arr, columns=["GROUP SIZE"])
        n2 = pd.DataFrame(n2_arr, columns=["GROUP SIZE"])
        n3 = pd.DataFrame(n3_arr, columns=["GROUP SIZE"])

        means1_arr = pd.concat([gene_of_study, y1, z1, n1], axis=1)
        means2_arr = pd.concat([y2, z2, n2], axis=1)
        means3_arr = pd.concat([y3, z3, n3], axis=1)

        all = pd.concat([means1_arr, means2_arr, means3_arr], axis=1)

        all_genes.append(gene_of_study)

        #     all_genes.append(genes)

        means1_table.append(means1_arr)
        means2_table.append(means2_arr)
        means2_table.append(means3_arr)
        # calculate the means std_dev and columns
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        # print(es_list)
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)

    return expressions_team1, expressions_team2, expressions_team3, means1_table, means2_table, means3_table


# print(expressions_team1)
# With this function we can conduct a meta-analysis (Random models, IV-Heg,and SMD)
# WARNING! : First we have to load our data (with the load_mage_data function) in order to execute this function
def calc_meta_data(expressions_team1, expressions_team2, expressions_team3, num_of_studies):
    global study, gene_list
    study = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
             'V', 'W',
             'X', 'Y', 'Z']
    all_list2 = []
    gene_list = []

    for i in range(len(expressions_team1)):
        #

        genes = expressions_team1[i][0]

        x1 = expressions_team1[i].iloc[:, 1:].mean(axis=1)
        x2 = expressions_team1[i].iloc[:, 1:].std(axis=1)

        n1 = len(ind1[i])
        n1 = [n1] * len(x1)

        n1 = pd.DataFrame(n1, columns=["GROUP SIZE"])

        x1 = pd.DataFrame(x1, columns=['MEANS ΤΕΑΜ 1'])
        x2 = pd.DataFrame(x2, columns=['STANDARD DEVIATIONS ΤΕΑΜ 1'])

        all_1 = pd.concat([genes, x1, x2, n1], axis=1)

        y1 = expressions_team2[i].iloc[:, 1:].mean(axis=1)
        y2 = expressions_team2[i].iloc[:, 1:].std(axis=1)

        n2 = len(ind2[i])
        n2 = [n2] * len(y1)

        z1 = expressions_team3[i].iloc[:, 1:].mean(axis=1)
        z2 = expressions_team3[i].iloc[:, 1:].std(axis=1)

        n3 = len(ind3[i])
        n3 = [n3] * len(z1)

        n2 = pd.DataFrame(n2, columns=["GROUP SIZE"])

        n3 = pd.DataFrame(n3, columns=["GROUP SIZE"])

        y1 = pd.DataFrame(y1, columns=['MEANS ΤΕΑΜ 2'])
        y2 = pd.DataFrame(y2, columns=['STANDARD DEVIATIONS ΤΕΑΜ 2'])

        z1 = pd.DataFrame(z1, columns=['MEANS ΤΕΑΜ 3'])
        z2 = pd.DataFrame(z2, columns=['STANDARD DEVIATIONS ΤΕΑΜ 3'])

        all_2 = pd.concat([y1, y2, n2], axis=1)

        all_3 = pd.concat([z1, z2, n3], axis=1)

        all_d = pd.concat([all_1, all_2, all_3], axis=1)
        all_list2.append(all_d)
    all_info = []
    all_genes2 = []
    for i in all_list2:
        all_info.append(i.values.tolist())
    for i in all_genes:
        all_genes2.append(i.values.tolist())

    items = []

    for j in all_info:
        items.append(j)
    items = list(np.concatenate(
        items))  # now everything is in the desired dataFrame,but we need to make the genes as keys in hash tables
    items = np.array(items)

    hash_table = []
    for item in items:
        hash_table.append({item[0]: item[1:]})
    hash_table = np.array(hash_table)
    # hash_table=list(hash_table.flat)

    res = defaultdict(list)
    for sub in hash_table:
        for key in sub:
            res[key].append(sub[key])
    hash_table = dict(res)
    res = {k: v for k, v in hash_table.items() if len(v) >= 2}
    hash_table = res

    # print(hash_table)

    x = pd.DataFrame(hash_table.values(), index=hash_table.keys()).T
    #     print(x)
    gene_list = list(x.columns)

    mult_var_list = []

    for gene in gene_list:
        #         print(gene)
        if len(x[gene]) < 2:
            gene_list.remove(gene)
    for gene in gene_list:
        for i in range(len(x[gene].dropna())):
            mult_args = np.array(x[gene][i]).astype(np.float)
            m1 = float(mult_args[0])
            m2 = float(mult_args[3])
            m3 = float(mult_args[6])

            st1 = float(mult_args[1])
            st2 = float(mult_args[4])
            st3 = float(mult_args[7])

            n1 = (mult_args[2])
            n2 = (mult_args[5])
            n3 = (mult_args[8])

            N = n1 + n2 + n3
            df = N - 3

            J = (math.gamma(df / 2) / (math.sqrt(df / 2) * math.gamma((df - 1) / 2)))

            Sp = math.sqrt(((n2 - 1) * st2 * st2 + (n1 - 1) * st1 * st1 + (n3 - 1) * st3 * st3) / (N - 3))

            d1 = (m2 - m1) / Sp
            d2 = (m3 - m1) / Sp

            var_d2 = (1 / n1) + (1 / n2) + (d1 ** 2) / (2 * N)
            var_d3 = (1 / n3) + (1 / n1) + (d2 ** 2) / (2 * N)

            g1 = J * d1
            g2 = J * d2

            var_g1 = (J ** 2) * var_d2
            var_g2 = (J ** 2) * var_d3

            cov_d1_d2 = (1 / n1) + (d1 * d2) / (2 * N)
            cov_g1_g2 = J * J * cov_d1_d2
            new_row = {"Gene": gene, "d1": d1, "d2": d2, "var_d2": var_d2, "var_d3": var_d3, "cov_d1_d2": cov_d1_d2,
                       "g1": g1, "g2": g2, "var_g1": var_g1, "var_g2": var_g2, "cov_g1_g2": cov_g1_g2, 'n1': n1,
                       'n2': n2, 'n3': n3}
            mult_var_list.append(new_row)

    df = pd.DataFrame(mult_var_list).T

    new_header = df.iloc[0]
    df = df[1:]  # take the data less the header row
    df.columns = new_header  # set the header row as the df header
    df = df.T

    new_df_list = []
    # unique genes and measure as a list
    for gene in gene_list:
        g1_list = list(df.loc[gene, 'g1'])
        g2_list = list(df.loc[gene, 'g2'])
        cov_g1g2_list = list(df.loc[gene, 'cov_g1_g2'])
        var_g2_list = list(df.loc[gene, 'var_g2'])
        var_g1_list = list(df.loc[gene, 'var_g1'])

        n1_list = list(df.loc[gene, 'n1'])
        n2_list = list(df.loc[gene, 'n2'])
        n3_list = list(df.loc[gene, 'n3'])

        w_list = []
        w2_list = []
        se_list = []
        es_list = []
        z_list = []
        w_i_list = []
        w_i_sqrt_list = []
        diff_list = []
        std_err_diff_list = []
        var_global1_list = []
        p_w_list = []

        for i in range(len(g1_list)):
            g1 = g1_list[i]
            g2 = g2_list[i]
            covg1g2 = cov_g1g2_list[i]
            varg2 = var_g2_list[i]
            varg1 = var_g1_list[i]

            n1 = n1_list[i]
            n2 = n2_list[i]
            n3 = n3_list[i]
            N = n1 + n2 + n3

            #             the W formula simplified in order to get the p_value  easily

            w = ((g1 ** 2) * varg2 + (- 2 * g1 * g2 * covg1g2) + (g2 ** 2) * varg1) / ((varg1 * varg2) - (covg1g2 ** 2))

            # me diff kai stdrer_diff kanw metanalisi   (global 2 - nadia )
            diff = g1 - g2  # TO thewrw effect size
            std_err_diff = math.sqrt(varg1 + varg2 - (2 * covg1g2))

            diff_list.append(diff)
            std_err_diff_list.append(std_err_diff)

            w2 = diff / std_err_diff
            w2_list.append(w2)
            # var_global2 = 2 * (1 - norm.cdf(abs(w2)))

            w_list.append(w)

            p_w = math.exp(-(w / 2))
            p_w_list.append(p_w)

            #             z-score
            #              z = −0.862 + √[0.743 − 2.404×log(P)]

            #             z = -0.862 + math.sqrt(0.743 - 2.404*math.log(p_w))

            # Monoplevros elegxos
            z = 1 - st.norm.ppf(p_w)
            z_list.append(z)
            # psakse ton diplebro

            w_i = math.sqrt(N)
            w_i_list.append(w_i ** 2)
            w_i_sqrt_list.append(math.sqrt(w_i) * z)

            se = math.sqrt((1 / n1) + (1 / n2) + (1 / n3))
            se_list.append(se)
            e_size = z * se
            es_list.append(e_size)

            var_global2 = 2 * (1 - norm.pdf(abs(w2)))

            var_global1_list.append(var_global2)

        # g1,g2 me ta standard erros tous  gia input
        # kai meta afou bgaloyme w kai d kanoyme metanalisi

        # Stoufer
        stoufer_weighted = sum(w_i_sqrt_list) / (math.sqrt(sum(w_i_list)))  # weighted stoufer
        stoufer = sum(z_list) / math.sqrt(len(z_list))  # stoufer

        p_stoufer_weighted = 2 * (1 - norm.cdf(abs(stoufer_weighted)))
        p_stoufer = 2 * (1 - norm.cdf(abs(stoufer)))

        # Fihser
        U = -2 * sum(np.log(p_w_list))
        fisher = U * sum([(- np.log(U) ** i) / math.factorial(i) for i in range(num_of_studies - 1)])
        p_fisher = 2 * (1 - chi2.cdf(abs(fisher), 2))
        # fisher_list.append(fisher)

        # Edgington1

        p_edg1 = (sum(p_w_list) ** num_of_studies) / math.factorial(num_of_studies)

        # Edgington2

        p_hat = sum(p_w_list) / num_of_studies

        U_edg2 = (0.5 - p_hat) * math.sqrt(12)
        p_edg2 = 2 * (1 - norm.cdf(abs(U_edg2)))

        new_row = {'Genes': gene, "g1": g1_list, "g2": g2_list, "cov_g1_g2": cov_g1g2_list, "varg1": var_g1_list,
                   'varg2': var_g2_list, 'w': w2_list, 'var_global_list': var_global1_list, "es": es_list,
                   'se': se_list,
                   'diff': diff_list, 'std_err_diff': std_err_diff_list,
                   'global1_p_fisher': p_fisher, 'global1_p_edg1': p_edg1, 'global1_p_edg2': p_edg2,
                   'global1_stoufer': p_stoufer,
                   'global1_stoufer_weighted': p_stoufer_weighted}

        new_df_list.append(new_row)
    new_df = pd.DataFrame(new_df_list)

    # print(new_df)

    final_df_list1 = []
    final_df_list2 = []
    final_df_list3 = []
    final_df_list4 = []

    for i in range(len(new_df.index)):
        es_g1 = list(new_df.loc[i, 'g1'])
        se_g1 = list(new_df.loc[i, "varg1"])
        es_g2 = list(new_df.loc[i, 'g2'])
        se_g2 = list(new_df.loc[i, "varg2"])

        es_g3 = list(new_df.loc[i, 'diff'])
        se_g3 = list(new_df.loc[i, "std_err_diff"])

        es_g4 = list(new_df.loc[i, 'es'])
        se_g4 = list(new_df.loc[i, "se"])

        # p_stoufer2= list(new_df.loc[i, "stoufer"])

        std_err_g1, mu_bar_g1, p_g1 = altmeta(es_g1, se_g1)
        std_err_g2, mu_bar_g2, p_g2 = altmeta(es_g2, se_g2)
        std_err_global2, mu_bar_global2, p_global2 = altmeta(es_g3, se_g3)
        std_err_global1, mu_bar_global1, p_global1 = altmeta(es_g4, se_g4)

        #         print(new_df['Gene'][i], Q, I2, tau2_DL, p_Q, se, z, mu_bar, p)

        new_row1 = {'Genes': new_df['Gene'][i], "g1": mu_bar_g1, 'se_g1': std_err_g1, 'p_g1': p_g1}
        final_df_list1.append(new_row1)

        new_row2 = {"g2": mu_bar_g2, 'se_g2': std_err_g2, 'p_g2': p_g2}
        final_df_list2.append(new_row2)

        new_row3 = {'global2': p_global2}
        final_df_list3.append(new_row3)

        new_row4 = {'global1_RE': p_global1}
        final_df_list4.append(new_row4)

    g1_data = pd.DataFrame(final_df_list1)
    g2_data = pd.DataFrame(final_df_list2)
    g3_data = pd.DataFrame(final_df_list3)
    g4_data = pd.DataFrame(final_df_list4)
    stoufer_data = new_df['global1_stoufer']
    stoufer_w_data = new_df['global1_stoufer_weighted']
    fisher_data = new_df['global1_p_fisher']
    edg1_data = new_df['global1_p_edg1']
    edg2_data = new_df['global1_p_edg2']

    final_df = pd.concat(
        [g1_data, g2_data, g4_data, stoufer_data, stoufer_w_data, fisher_data, edg1_data, edg2_data, g3_data], axis=1)

    # print(final_df.head())
    return final_df


def run(settings, data, filepath):
    alpha = float(settings['alpha'])
    venn_correction = settings['venn_correction']
    venn_choice = settings['venn_choice']
    multiple_tests = settings['multiple_tests']
    num_of_studies = len(data)
    expressions_team1, expressions_team2, expressions_team3, means1_table, means2_table, means3_table = split_data(data)
    df = calc_meta_data(expressions_team1, expressions_team2, expressions_team3, num_of_studies)
    df = df.dropna()  # remove Nan lines

    if (multiple_tests != 'none') & (venn_correction != 'none'):
        if multiple_tests == 'one_step':
            one_df1 = get_one_step_methods(df, alpha, 'p_g1')
            one_df2 = get_one_step_methods(df, alpha, 'p_g2')
            one_df3 = get_one_step_methods(df, alpha, venn_choice)
            one_df = pd.concat([one_df1, one_df2, one_df3], axis=1)
            df = (pd.concat([df, one_df], axis=1))

        if multiple_tests == 'step_down':
            one_df1 = get_step_down_methods(df, alpha, 'p_g1')
            one_df2 = get_step_down_methods(df, alpha, 'p_g2')
            one_df3 = get_step_down_methods(df, alpha, venn_choice)
            one_df = pd.concat([one_df1, one_df2, one_df3], axis=1)
            df = (pd.concat([df, one_df], axis=1))

        if multiple_tests == 'step_up':
            one_df1 = get_step_up_methods(df, alpha, 'p_g1')
            one_df2 = get_step_up_methods(df, alpha, 'p_g2')
            one_df3 = get_step_up_methods(df, alpha, venn_choice)
            one_df = pd.concat([one_df1, one_df2, one_df3], axis=1)
            df = (pd.concat([df, one_df], axis=1))

        df.to_csv("results_multivariate.txt", sep='\t', mode='w')
        df = df.loc[:, ~df.columns.duplicated()]
        genes_venn = list(df['Genes_' + multiple_tests])

        list1 = list(df[venn_correction + '_p_g1'])
        list2 = list(df[venn_correction + '_p_g2'])
        list3 = list(df[venn_correction + "_" + venn_choice])

        l1 = []
        l2 = []
        l3 = []

        for i in range(len(genes_venn)):
            if list1[i] == 1:
                l1.append(genes_venn[i])
            if list2[i] == 1:
                l2.append(genes_venn[i])

            if list3[i] == 1:
                l3.append(genes_venn[i])

    genes_venn = list(df['Genes'])

    list1 = list(df['p_g1'])
    list2 = list(df['p_g2'])
    list3 = list(df[venn_choice])

    l1 = []
    l2 = []
    l3 = []

    for i in range(len(genes_venn)):
        if list1[i] < 0.05:
            l1.append(genes_venn[i])
        if list2[i] < 0.05:
            l2.append(genes_venn[i])
        if list3[i] < 0.05:
            l3.append(genes_venn[i])

    if settings.get('plots') == 'YES':
        plots.multivariate_plots(l1, l2, l3, venn_correction, venn_choice, filepath)

    return df
