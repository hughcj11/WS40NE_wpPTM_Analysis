# PTM_Proteomic_Analysis
https://github.com/hughcj11/WS40NE_wpPTM_Analysis

The following code works on a csv file downloaded from Skyline, which should contain the following information in this order: Protein Name, Protein Accession, Protein Sequence Coverage,Protein Sequence, Protein Description, and the Normalized Area of each Replicate. If Skyline is listing replicates one by one in a row instead of columns, check the "Pivot Replicates" button. It is critical that the csv file is ordered this way so that the code reads the correct columns.
  
To run your specific file, you must first update the base path so that the code can locate your specific data file. Once the base path is updated, all files generated will appear in the same folder as your Skyline data.


Protein Abundance Analysis: (code files: Cell Proteomics_Abundance)
 The following code will generate the following documents using the Skyline csv file:

     quantile: quantile normalization of data

     SighProteinAbundance: lists all the sig proteins 

     SighProteinAbundanceEnsembleID:lists all the ensemble ids for the sig proteins 

     ProteinsforCLUST: Lists all proteins with their matching EnsembleID.


R Code for Figures:
The code "Cell Proteomics R Code" is used in R to generate abundance maps/PCA plots. 
