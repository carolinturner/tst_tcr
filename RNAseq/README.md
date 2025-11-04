# Scripts for RNAseq analyses

The scripts provided here were used for RNAseq analyses in the TST TCR manuscript.
_RNAseq_scripts_ generate the data required for _Plotting_scripts_ to make Figure 1 and Supplementary Figure 2. _RNAseq_Table_script_ makes Table 1.

At time of peer-reviewed publication of the manuscript, processed RNAseq source data can be downloaded from [Array Express](https://www.ebi.ac.uk/arrayexpress), accession number E-MTAB-14687.
 
Metadata = `sdrf.tsv`
Raw counts matrix = `rawcounts.csv`
TPM matrix = `tpm_PC0.001_log2_genesymbol_dedup.csv`

Overview of analysis scripts:
* `RNAseq_script_1.R`: Prepare data and differential gene expression analyses
* `RNAseq_script_2.R`: SARtools DeSeq2 analysis (TST_D2 vs saline) and prep for subsequent analysis
* `RNAseq_script_3.R`: SARtools DeSeq2 analysis (TST_D7 vs saline) and prep for subsequent analysis
* `RNAseq_script_4.R`: SARtools DeSeq2 analysis (all TST vs saline) and prep for subsequent analysis
* `RNAseq_script_5.R`: SARtools DeSeq2 analysis (D7 vs D2 within integrated TST response) and prep for subsequent analysis
* `RNAseq_script_6.R`: XGR pathway analysis
* `RNAseq_script_7.R`: Correlation analysis of upstream regulator target genes (TST_D2 vs saline) and prep for Gephi visualisation
* `RNAseq_script_8.R`: Correlation analysis of upstream regulator target genes (TST_D7 vs saline) and prep for Gephi visualisation
* `RNAseq_script_9.R`: Correlation analysis of upstream regulator target genes (all TST vs saline) and calculation of module scores per sample
* `RNAseq_script_10.R`: Module analyses

To make:
* **Table 1**: RNAseq_Table_script.R
* **Figure S2A**: RNAseq_script_1 &rarr; RNAseq_script_2 &rarr; IPA upstream regulator analysis &rarr; RNAseq_script_7 &rarr; Gephi network visualisation
* **Figure S2B**: RNAseq_script_1 &rarr; RNAseq_script_3 &rarr; IPA upstream regulator analysis &rarr; RNAseq_script_8 &rarr; Gephi network visualisation
* **Figure S2C-D**: RNAseq_script_1 &rarr; RNAseq_script_4 &rarr; IPA upstream regulator analysis &rarr; RNAseq_script_9 &rarr; Plotting_script_FigureS2C-D
* **Figure S2E**: RNAseq_script_1 &rarr; RNAseq_script_5 &rarr; Plotting_script_FigureS2E
* **Figure S2F**: RNAseq_script_1 &rarr; RNAseq_script_2 and RNAseq_script_3 &rarr; RNAseq_script_5 &rarr; RNAseq_script_6 &rarr; Plotting_script_FigureS2F
* **Figure 1A-B**: RNAseq_script_1 &rarr; RNAseq_script_10 &rarr; Plotting_script_Figure1A-B
* **Figure 1C**: RNAseq_script_1 &rarr; RNAseq_script_10 &rarr; Plotting_script_Figure1C