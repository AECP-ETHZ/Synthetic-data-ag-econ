# Synthetic-data-ag-econ

This repository provides the R codes for the publication "A note on synthetic data for replication purposes in agricultural economics".

Please use the following citation: Wimmer, S., Finger, R. (2023). "A note on synthetic data for replication purposes in agricultural economics." Journal of Agricultural Economics. 

The files are as follows (note that only 01_Analysis_Norway.R can be readily run without further data access):

- 01_Analysis_Norway.R: The program stores all tables and figures based on the Norwegian dairy production data presented in the article in the subfolders 'Tables/Norway' and 'Figures/Norway'. The data are automatically downloaded from the internet. 

- 02_Prepare_FADN.R: This program prepares the FADN data for production function and frontier estimation. The data (stored in the subfolder 'Data') are not publicly available but can be requested from DG-Agri. 

- 03_Analysis_FADN.R: The program stores all tables and figures based on the FADN data for German crop farms presented in the appendix in the subfolders 'Tables/FADN' and 'Figures/FADN'. The data are not publicly available but can be requested from DG-Agri.

The script was successfully run with R version 4.2.0.