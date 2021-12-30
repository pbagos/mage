import numpy as np
import scipy.stats as stats
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn3_circles, venn3


def meta_analysis_plots(meta_analysis_df, filepath):
    alpha = 0.01

    # QQ plot
    plt.figure("QQ plot ")
    measurements = meta_analysis_df["Effect size (Hedge's g)"]
    tau = meta_analysis_df['Tau_Squared']
    q = meta_analysis_df['Q']
    i2 = meta_analysis_df['I_Squared']
    stats.probplot(measurements, dist="norm", plot=plt)
    plt.savefig(filepath + '/qq_plot.png')

    # Tau2 plot
    plt.figure("Tau - Squared DL Histogram ")
    plt.hist(tau, color='r', density=True, edgecolor='black', linewidth=1.2)
    plt.title("Tau - Squared DL Histogram")
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.savefig(filepath + '/tau2.png')

    # Q plot
    plt.figure("Q Histogram ")
    plt.hist(q, color='g', density=True, edgecolor='black', linewidth=1.2)
    plt.title("Q Histogram")
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.savefig(filepath + '/q.png')

    # I2 plot
    plt.figure("I - Squared Histogram ")
    plt.hist(i2, color='b', density=True, edgecolor='black', linewidth=1.2)
    plt.title("I - Squared DL Histogram ")
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.savefig(filepath + '/i2.png')

    # Volcano Plot
    df = meta_analysis_df[meta_analysis_df.p_value != 0]
    p_values = np.array(df['p_value'].to_numpy().tolist(), dtype=float)

    smd = df["Effect size (Hedge's g)"]
    log_p = -(np.log10(p_values))

    color = []

    for i in range(len(p_values)):
        if p_values[i] < alpha:
            color.append(2)
        else:
            color.append(1)

    plt.figure("Volcano plot")
    plt.scatter(smd, log_p, c=color, cmap='RdYlGn')
    plt.title("Volcano Plot")
    plt.xlabel("Effect Size Hedge's g")
    plt.ylabel("- log10 p_values")
    plt.savefig(filepath + "/volc_plot.png")

    # plt.show()


def multivariate_plots(list1, list2, list3, venn_correction, venn_choice, filepath):
    set1 = set(list1)
    set2 = set(list2)
    set3 = set(list3)
    plt.title("p-values Venn diagram" + " - " + venn_correction + " correction")
    venn3([set1, set2, set3], ('g1', 'g2', venn_choice))
    venn3_circles(subsets=[set1, set2, set3], linestyle='dashed')
    plt.savefig(filepath + "/venn_plot.png")


def ea_manhattan_plot(data, filepath, p_threshold):
    df = data.loc[data.p_value < p_threshold, ['source', 'native', 'p_value', 'negative_log10_of_adjusted_p_value']]
    df = df.where(df ['negative_log10_of_adjusted_p_value'] <= 16)
    # Generate Manhattan plot: (#optional tweaks for relplot: linewidth=0, s=9)
    df['group_id'] = df.groupby('source').ngroup()
    df.sort_values(['group_id', 'native'], inplace=True)
    df.reset_index(inplace=True, drop=True)
    df['i'] = df.index

    plt.subplots(dpi=300)
    plt.figure("Manhattan plot")

    # palette: Paired, Set2, bright, dark, coloblind, deep
    plot = sns.relplot(data=df, x='i', y='negative_log10_of_adjusted_p_value', aspect=4,
                       hue='source', palette='bright', legend='full')
    plot.ax.axhline(16, linestyle='--', linewidth=1)

    plot.set_axis_labels('chromosomes', '-log10(Padj)')
    plot.fig.suptitle('Enrichment Analysis - Manhattan plot', fontsize=16)
    plot.ax.set_xticks(df.groupby('group_id')['i'].median())
    plot.ax.set_xticklabels(df['source'].unique(), rotation=25)
    plt.savefig(filepath + "/manhattan_plot.png")


def ea_heatmap_plot(data, filepath):
    # Get unique sources
    sourses = data.source.unique()

    # Create a heatmap plot for each source
    for source in sourses:
        genes_all = []
        for index, row in data[data.source == source].iterrows():
            genes_all.extend(row['intersections'])

        genes_all = list(set(genes_all))
        # print(genes_all)
        df = pd.DataFrame(columns=genes_all)

        for index, row in data[data.source == source].iterrows():
            df.loc[row['native'], row['intersections']] = row['negative_log10_of_adjusted_p_value']

        df = df.fillna(0)
        # print(df)
        plt.figure("Heatmap plot - " + source, figsize=(30, 15), dpi=80)
        sns.heatmap(df, cmap='RdYlGn_r')
        sns.set(font_scale=1.2)

        plt.title(source, fontsize=20)
        plt.xlabel('Genes', fontsize=15)  # x-axis label with fontsize 15
        plt.ylabel('Terms', fontsize=15)  # y-axis label with fontsize 15
        plt.savefig(filepath + "/heatmap_plot_" + source + ".png")
        # plt.show()