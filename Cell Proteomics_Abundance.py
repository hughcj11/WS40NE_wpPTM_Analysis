#Author: Chelsea Hughes
import numpy as np
from pandas import read_csv
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from numpy import inf

#Since the log2 did not normalize data, quantile normalization was used
data = read_csv("/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Protein Abundance/Protein Abundance for WS40NEwpPTM.csv",header=0)
datastats2= data.iloc[:,[0,4]]
datastats= data.iloc[:,5:51]
def quantileNormalize(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df
result = quantileNormalize(datastats)
result2 = pd.concat([datastats2, result], axis=1)
result2.to_csv('/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Protein Abundance/quantile.csv', index=False) 

####T-test on normalized relative abundance
#################Normoxic vs Recovery

data = read_csv('/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Protein Abundance/quantile.csv',header=0)
datastats= data.iloc[:,0:2]
datastats["Normoxic_Average"]=(data.iloc[:,14:26].sum(axis=1)/12)*100
datastats["Recovery"]=(data.iloc[:,26:38].sum(axis=1)/12)*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['Recovery'] / (datastats['Normoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data.iloc[:,14:26], data.iloc[:,26:38],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
# Write DataFrame to CSV
datastats.to_csv('/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Protein Abundance/DatastatsNormoxicVRecovery.csv', index=False)


#################Normoxic vs 4 days anoxia

data = read_csv('/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Protein Abundance/quantile.csv',header=0)
datastats= data.iloc[:,0:2]
datastats["Normoxic_Average"]=(data.iloc[:,14:26].sum(axis=1)/12)*100
datastats["4dAnoxia"]=(data.iloc[:,2:14].sum(axis=1)/12)*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['4dAnoxia'] / (datastats['Normoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data.iloc[:,14:26], data.iloc[:,2:14],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
# Write DataFrame to CSV
datastats.to_csv('/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Protein Abundance/DatastatsNormoxicV4dAnoxia.csv', index=False)


#################4Danoxia vs Recovery

data = read_csv('/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Protein Abundance/quantile.csv',header=0)
datastats= data.iloc[:,0:2]
datastats["Recovery"]=(data.iloc[:,26:38].sum(axis=1)/12)*100
datastats["4dAnoxia"]=(data.iloc[:,2:14].sum(axis=1)/12)*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['Recovery'] / (datastats['4dAnoxia']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data.iloc[:,2:14], data.iloc[:,26:38],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
# Write DataFrame to CSV
datastats.to_csv('/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Protein Abundance/DatastatsRecoveryV4DAnoxia.csv', index=False)



#Code to find all sig proteins and compile them in one document

DataRecv4d= read_csv('/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Protein Abundance/DatastatsRecoveryV4DAnoxia.csv') #
DataNormv4d= read_csv('/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Protein Abundance/DatastatsNormoxicV4dAnoxia.csv')#
DataNormvRec= read_csv('/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Protein Abundance/DatastatsNormoxicVRecovery.csv')

DataNormv4dFiltered= DataNormv4d[DataNormv4d['corrected_p_values']<.05]
DataNormv4dFiltered.rename(columns={"log2FC": "Normoxia vs 4d Anoxia log2FC"}, inplace=True) 
DataNormv4dFiltered= DataNormv4dFiltered.drop(columns=["Normoxic_Average", "4dAnoxia",'corrected_p_values','T-test_t_statistic', 'T-test_p_value'])

DataNormvRecFiltered= DataNormvRec[DataNormvRec['corrected_p_values']<.05]
DataNormvRecFiltered.rename(columns={"log2FC": "Normoxia vs Recovery log2FC"}, inplace=True) 
DataNormvRecFiltered= DataNormvRecFiltered.drop(columns=["Normoxic_Average","Recovery",'corrected_p_values','T-test_t_statistic', 'T-test_p_value'])

DataRecv4dFiltered= DataRecv4d[DataRecv4d['corrected_p_values']<.05]
DataRecv4dFiltered.rename(columns={"log2FC": "4d Anoxia vs Recovery log2FC"},inplace=True) 
DataRecv4dFiltered= DataRecv4dFiltered.drop(columns=["4dAnoxia","Recovery",'corrected_p_values','T-test_t_statistic', 'T-test_p_value'])
merged = pd.merge(DataRecv4dFiltered, DataNormv4dFiltered, on=["Protein Name", "Protein Description"], how='outer')
merged = pd.merge(merged, DataNormvRecFiltered, on=["Protein Name", "Protein Description"], how='outer')

merged =merged[["Protein Name", "Protein Description","Normoxia vs 4d Anoxia log2FC", "Normoxia vs Recovery log2FC", "4d Anoxia vs Recovery log2FC"]]


merged.to_csv('/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Protein Abundance/SighProteinAbundance.csv', index=False)




#Code to find all the ensemble ids for the sig proteins 
SigProteins = pd.read_csv("/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Protein Abundance/SighProteinAbundance.csv")

EnsembleID = pd.read_csv("/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Testing code/Orthofinder_blast_merge_Hs_230628_MapIDS.csv", encoding='latin-1')

# Merge the two dataframes based on the "Protein Name" column
merged_df = pd.merge(SigProteins, EnsembleID, how='left', left_on='Protein Name', right_on='Austrofundulus_limnaeus')

# Add a new column to the first dataframe with values from the "EnsemblID_Hs" column
merged_df['EnsemblID_Hs'] = merged_df['EnsemblID_Hs']

# Write the merged dataframe to a new CSV file
merged_df.to_csv("/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Protein Abundance/SighProteinAbundanceEnsembleID.csv")




#Code to find all the ensemble ids for all proteins to use in CLUST
SigProteins = pd.read_csv('/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Protein Abundance/quantile.csv')

EnsembleID = pd.read_csv("/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Testing code/Orthofinder_blast_merge_Hs_230628_MapIDS.csv", encoding='latin-1')

# Merge the two dataframes based on the "Protein Name" column
merged_df = pd.merge(SigProteins, EnsembleID, how='left', left_on='Protein Name', right_on='Austrofundulus_limnaeus')

# Add a new column to the first dataframe with values from the "EnsemblID_Hs" column
merged_df['EnsemblID_Hs'] = merged_df['EnsemblID_Hs']

# Write the merged dataframe to a new CSV file
merged_df.to_csv("/Users/chelseahughes/Desktop/PTM Proteomic Analysis/Protein Abundance/ProteinsforCLUST.csv")