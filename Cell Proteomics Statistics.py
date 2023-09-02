#Author: Chelsea Hughes
import csv
import numpy as np
from pandas import read_csv
import pandas as pd
from scipy.stats import f_oneway
from scipy.stats import ttest_rel
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from numpy import inf
base_path = "/Users/chelseahughes/Desktop/Histone Analysis/Calculation for Embryo Samples/"

#####Generating values for analysis
##Unique hPTMs
#The below code generates the relative abundance value in a csv file
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/Embryo Library/Replicate calculations.csv",header=0)
datastats= data.iloc[:,0:8]
sample=[f"Sample {i+1}" for i in range(0,12)]
datastats[sample]= (data.iloc[:,9:21].to_numpy() /(data.iloc[:,22:34])*100).to_numpy()
datastats=datastats.replace(np.nan, 0)
datastats.to_csv(base_path+'/Full_Relative_Abundance.csv', index=False)
#The below code generates the beta value in a csv file
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/Embryo Library/Replicate calculations.csv",header=0)
datastats= data.iloc[:,0:8]
sample=[f"Sample {i+1}" for i in range(0,12)]
datastats[sample]= data.iloc[:,9:21].to_numpy() /(data.iloc[:,22:34]+100).to_numpy()
datastats=datastats.replace(-np.inf, 0)
datastats.to_csv(base_path+'/Full_beta_values.csv', index=False)


#The below code generates the m value in a csv file
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/Embryo Library/Replicate calculations.csv",header=0)
datastats= data.iloc[:,0:8]
sample=[f"Sample {i+1}" for i in range(0,12)]
datastats[sample]=np.log2((data.iloc[:,9:21].to_numpy()/((data.iloc[:,22:34]+100).to_numpy()))/(1-(data.iloc[:,9:21].to_numpy()/((data.iloc[:,22:34]+100).to_numpy()))))
datastats=datastats.replace(-np.inf, 0)
datastats.to_csv(base_path+'/Full_M_Values.csv', index=False)

##Global calculations
#The below code generates the global relative abundance value in a csv file
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/Embryo Library/Replicate_calculations_Total_PTMs.csv",header=0)
datastats= data.iloc[:,0:1]
sample=[f"Sample {i+1}" for i in range(0,12)]
datastats[sample]= (data.iloc[:,1:13].to_numpy()/(data.iloc[:,14:26])*100).to_numpy()
datastats=datastats.replace(-np.inf, 0)
datastats.to_csv(base_path+'/Embryo Anoxia Total PTM.csv', index=False)


####Statistical analysis
##Unique hPTMs
#The below code performs statistics on the relative abundance data for each unique hPTM (residue,PTM,and histone)
data2 = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculation for Embryo Samples/Full_Relative_Abundance.csv",header=0)
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/Embryo Library/Replicate calculations.csv",header=0)
datastats= data.iloc[:,0:8]
datastats["Normoxic_Average"]=(data.iloc[:,9:15].sum(axis=1)/data.iloc[:,22:28].sum(axis=1))*100
datastats["Anoxic_Average"]=(data.iloc[:,15:21].sum(axis=1)/data.iloc[:,28:34].sum(axis=1))*100
datastats=datastats.replace(np.inf, 0)
datastats['log2FC'] = np.log2(datastats['Anoxic_Average'] / (datastats['Normoxic_Average']+.0000001)) #.0000001 added to remove zeroes
datastats=datastats.replace(-np.inf, 0)
# datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:14],data.iloc[:,14:20],axis=1)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,8:14], data2.iloc[:,14:20],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]

# Write DataFrame to CSV
datastats.to_csv(base_path+'/Datastats.csv', index=False)

##Global statistics
#The below code performs statistics on the relative abundance data global hPTM changes (independent of residue and histone)
data2 = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/Calculation for Embryo Samples/Embryo Anoxia Total PTM.csv",header=0)
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/Embryo Library/Replicate_calculations_Total_PTMs.csv",header=0)
datastats= data.iloc[:,0:1]
datastats["Normoxic_Average"]=(data.iloc[:,1:7].sum(axis=1)/data.iloc[:,14:20].sum(axis=1))*100
datastats["Anoxic_Average"]=(data.iloc[:,7:13].sum(axis=1)/data.iloc[:,20:27].sum(axis=1))*100
datastats['log2FC'] = np.log2(datastats['Anoxic_Average'] / (datastats['Normoxic_Average']+.0000001)) #.0000001 added to remove zeroes
# datastats['Fvalue'],datastats['Anova_Pvalue'] = f_oneway(data.iloc[:,8:14],data.iloc[:,14:20],axis=1)
datastats['T-test_t_statistic'], datastats['T-test_p_value'] = ttest_ind(data2.iloc[:,1:7], data2.iloc[:,7:13],axis=1)
datastats['T-test_p_value'] = datastats['T-test_p_value'].fillna(1)
datastats['corrected_p_values'] = multipletests(datastats['T-test_p_value'], method='fdr_bh')[1]
datastats.to_csv(base_path+'/DatastatsGlobal.csv', index=False)

#The below document shows the relative coverage of each modifiable residue (a residue shown as capable of having a PTM). For example, how often K covered by PTMs?
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/Embryo Library/Replicate calculations.csv",header=0)
datastats= pd.DataFrame()
GroupedAA=data.groupby('Amino Acid').sum().reset_index()
datastats['Amino Acid']=GroupedAA['Amino Acid']
datastats["Normoxic_Average"]=(GroupedAA.iloc[:,9:15].sum(axis=1)/GroupedAA.iloc[:,22:28].sum(axis=1))*100
datastats["Anoxic_Average"]=(GroupedAA.iloc[:,15:21].sum(axis=1)/GroupedAA.iloc[:,28:34].sum(axis=1))*100
datastats.to_csv(base_path+'/ResidueCoverage.csv', index=False)

#The below document shows the relative coverage by a specific PTM for each modifiable residue (a residue shown as capable of having a PTM). For example, how often K covered by Ub?
data = read_csv("/Users/chelseahughes/Desktop/Histone Analysis/code/Embryo Library/Replicate calculations.csv",header=0)
datastats= pd.DataFrame()
GroupedAA=data.groupby(['Amino Acid',"PTM Description"]).sum().reset_index()
datastats['Amino Acid']=GroupedAA['Amino Acid']
datastats["PTM Description"]=GroupedAA["PTM Description"]
#Fix this- the math below divides only by where it could be based on where that modification shows up
datastats["Normoxic_Average"]=(GroupedAA.iloc[:,9:15].sum(axis=1)/GroupedAA.iloc[:,22:28].sum(axis=1))*100
datastats["Anoxic_Average"]=(GroupedAA.iloc[:,15:21].sum(axis=1)/GroupedAA.iloc[:,28:34].sum(axis=1))*100
datastats.to_csv(base_path+'/ResidueCoverageByPTM.csv', index=False)
  
  