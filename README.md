# LA_ICP_MS-Routines_FMNH
Routines to read and analyze LA-ICP-MS data from the Field Museum's Elemental Analysis Facility

Significant contributions include:
- Reading in large amounts of cleaned LA-ICP-MS data from xlsx files with each day as a separate sheet
- Tidying the LA-ICP-MS data into a tidy data frame
- Functions for converting elements calculated as % oxide to ppm
- Writing the file to a csv in preparation for data anlaysis

These routines and the analysis are very much in progress as of April 2018. 

Features coming include:

- Conversion of data to log scale
- Dendrogram object creation and analysis
- Plotting, coloring, and creating cutoffs for dendrograms as a heuristic tool
- Creating clusters based on dendrogram objects
- Principal Components Analysis using clusters from dendrogram objects
- Building scree plots to show variance explained
- Exploring eigenvector loadings on PCA
- Using Partitioning Around Mediods based on Silhouette Analysis to compare to HCA clusters
- Visualizing a prior information about sherds to provide insight to the analysis 

