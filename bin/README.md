# WHAT, WHO, WHEN

- Folder containing the files associated with Lee et al. Predicting the capsid architecture of tailed phages from metagenomic data

- Diana Y. Lee, Luque Lab, SDSU / dlee@sdsu.edu
- 2021-09-24, revised 2021-11-12


# FOLDERS & FILES

**FILE:** Figure_3_Statistical_Model_Updated.ipynb
--> Jupyter notebook for Python 3.0 
Creates a linear regression model for relating the genome length of a tailed phage to the capsid architecture 
as measured by the T-number.
Evaluates the accuracy of that model.
Creates the base illustrations for Figure 2 of the paper
--> Requires (data folder):
HRphagedata.csv
--> Creates (results folder):
    Fig3_DNA_vs_T-number(ticks).svg
    Fig3_DNA_vs_T-number.svg .png
    Fig3_DNA_vs_T_model_Rarification.svg .png
    Figure3_update.db (saved kernel state)

**FILE:** Figure_4_model_application_to_MCP_Genome_DB_Updated.ipynb
--> Jupyter notebook for Python 3.0 
Creates a kernel distribution for analyzing the MCP database
Applies the G2T model to the MCP database
--> Requires:
AllDNA.csv
phage_functions.ipynb
PHAGE_TABLE3.xlsx
--> Creates:
Fig4_MCP_Kernel_Density_all.svg
Fig4a_MCP_Kernel_Density.svg / .png
Fig4a_MCP_Kernel_Density_shaded.svg / .png
Fig4b_Predicted_T-numbers_for_MCP.svg / .png
Fig4c_Percent_architectures_predicted_by_T-number.svg / .png
#MCP_Kernel_Density_(hex).svg
#MCP_Kernel Density (trihex).svg
#Fig4bSI-Predicted T-numbers for Unassigned MCP.svg / .png
**Database file:** (saved kernel state)
Figure4_update.db

# HOW
