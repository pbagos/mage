import pandas as pd
import numpy as np



def run(settings, study, platform):
        transformation_method = settings['transformation_method']
        if settings['gene_data_online'] == 'YES':
            # Load data from web
            try:
                gene_history = pd.read_csv("https://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz", delimiter="\t",
                                           usecols=['GeneID', 'Discontinued_GeneID'])

                print(gene_history.head())
            except:
                print("Error load Gene History data from web")
                gene_history = pd.read_csv(settings['gene_history_file'], delimiter="\t")

            try:
                NCBI_intel = pd.read_csv("https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
                                         delimiter="\t", usecols=['GeneID', 'Symbol'])
            except:
                print("Error load data from web")
                NCBI_intel = pd.read_csv(settings['homo_sapiens_file'], delimiter="\t")

        else:
            # Load Data from local folder
            gene_history = pd.read_csv(settings['gene_history_file'], delimiter="\t")
            NCBI_intel = pd.read_csv(settings['homo_sapiens_file'], delimiter="\t")

        Platform_intel = pd.read_csv(settings['platforms_folder'] + platform + ".txt", delimiter="\t")

        # Get dataframe without first line
        data = study.iloc[2:].copy()
        data.rename(columns={data.columns[0]: "ID_REF"}, inplace=True)

        # Get Gene IDs
        dataIDs = pd.DataFrame(data["ID_REF"])

        # Gene name update based on gene_history
        probe_ID_platform = pd.merge(dataIDs, Platform_intel, how="inner", left_on="ID_REF", right_on="ID").drop(
            columns="ID")

        probe_ID_gene_history = pd.merge(probe_ID_platform, gene_history, how="inner", left_on="SPOT_ID",
                                         right_on="Discontinued_GeneID").drop(columns=["Discontinued_GeneID", "SPOT_ID"])
        probe_ID_gene_history_UPDATED = pd.merge(probe_ID_platform, probe_ID_gene_history, how="outer", left_on="ID_REF",
                                                 right_on="ID_REF")
        probe_ID_gene_history_UPDATED['SPOT_ID'] = np.where(probe_ID_gene_history_UPDATED.GeneID.notna(),
                                                            probe_ID_gene_history_UPDATED['GeneID'],
                                                            probe_ID_gene_history_UPDATED['SPOT_ID'])

        probe_ID_gene_history_UPDATED.drop(columns="GeneID", inplace=True)

        probe_ID_final = probe_ID_gene_history_UPDATED[probe_ID_gene_history_UPDATED['SPOT_ID'] != '-'].copy()
        probe_ID_final.dropna(axis=0, how="any", inplace=True)
        probe_ID_final['SPOT_ID'] = probe_ID_final['SPOT_ID'].astype(int)

        probe_ID_NCBI = pd.merge(probe_ID_final, NCBI_intel, how="inner", left_on="SPOT_ID", right_on="GeneID").drop(
            columns="GeneID")

        final_file = pd.merge(probe_ID_NCBI, data, how="inner", left_on="ID_REF", right_on="ID_REF")
        final_file.drop(columns="ID_REF", inplace=True)

        for col in final_file.columns[2:]:
            final_file[col] = final_file[col].astype(float)

        fun = {i: transformation_method for i in list(final_file.columns.values)[2:]}
        output = final_file.groupby(by=['SPOT_ID', 'Symbol']).agg(fun).reset_index()
        output.drop(columns="SPOT_ID", inplace=True)

        output.loc[-1] = study.iloc[1] # adding a row
        output.loc[-2] = study.iloc[0]  # adding a row
        output.index = output.index + 2  # shifting index
        output.sort_index(inplace=True)

        return output




def run_updated_genes(settings, study):
    transformation_method = settings['transformation_method']

    data = study.iloc[2:].copy()
    data.rename(columns={data.columns[0]: "ID_REF"}, inplace=True)
    # Get Gene IDs
    dataIDs = pd.DataFrame(data["ID_REF"])

    for col in data.columns[2:]:
        data[col] = data[col].astype(float)

    fun = {i: transformation_method for i in list(data.columns.values)[2:]}

    output = data.groupby(by=['ID_REF']).agg(fun).reset_index()

    # output['ID_REF'][0] = 'ID'
    # output['ID_REF'][1] = 'CLASS'
    output.loc[-1] = study.iloc[1] # adding a row
    output.loc[-2] = study.iloc[0]  # adding a row
    output.iloc[-1,0] = 'ID'
    output.iloc[-2,0] = 'CLASS'
    output.index = output.index + 2  # shifting index
    output.sort_index(inplace=True)
    output.columns = range(output.columns.size)
    print(output.head())
    return output
