# PCA
This repo is currently under construction.

# Make PC Loadings
- From reference data, extract sites in reference data that are common (MAF>0.02) present in TOPMed, and in approximate linkage equilibrium
- Remove related individuals (<3rd degree, using --king-cutoff in plink2) from refernece data
- Run PCA, save loadings
