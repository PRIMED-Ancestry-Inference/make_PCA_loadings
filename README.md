# PCA
This repo is currently under construction.

2 wdls available in repo:
- create PC loadings from reference data
- project samples onto PC loadings

# Create PC Loadings
- From reference data, extract sites in reference data that are common (MAF>0.02) present in TOPMed, and in approximate linkage equilibrium
- Remove related individuals (<3rd degree, using --king-cutoff in plink2) from refernece data
- Run PCA, save loadings

# Project Samplels on PC's
- Extract sites in PC loadings, require that __% of sites from loadings be present in data or fail
- Use plink2 --score
