# Scripts for TCRseq analyses

The scripts provided here were used for TCRseq analyses in the [TST TCR publication](link here once available). `TCRseq_scripts` generate the data required for `Plotting_scripts` to make Figures 2-5 and Supplementary Figures S3-S11 as detailed below. Figures 4 and S7 were assembled in Inkscape.

**Input data:**
* Processed TCRseq data and metadata can be downloaded from [UCL's Research Data Repository](add here link to DOI 10.5522/04/28049606 once public).
* Antigen-reactive CDR3 sequences collated from published sources are available as Table S2 of the [manuscript](link here once available).
* HLA imputations for participants are available as Table S3.
* Metaclone summaries (including HLA association and regular expression patterns) are available as Tables S4-S7. These metaclones were discovered using the code available in the metaclone directory of this GitHub repository.
* All source data should be saved in a subdirectory called `data`.
* Supplementary tables should be saved as individual .csv files, e.g. `TableS2.csv`

**Overview of scripts:**
* `TCRseq_script_1.R`: Pre-processing of data
* `TCRseq_script_2.R`: Calculate diversity indices
* `TCRseq_script_3.R`: Calculate abundance of published antigen-reactive TCRs in blood and TST
* `TCRseq_script_4.ipynb`: Coincidence analysis (down-sampled repertoires)
* `TCRseq_script_5.ipynb`: Coincidence analysis (full repertoires)
* `TCRseq_script_6.R`: Define private antigen-reactive TCRs
* `TCRseq_script_7.R`: Calculate abundance of private antigen-reactive TCRs in blood and TST
* `TCRseq_script_8.R`: Expansion of antigen-reactive TCRs in paired TST samples from day 2 and day 7
* `TCRseq_script_9.R`: Calculate abundance of metaclones in blood and TST 
* `TCRseq_script_10.R`: Identify most abundant and most public metaclones in day 7 TST
* `TCRseq_script_11.R`: Calculate abundance (and odds ratios) of metaclones and gliph2 patterns in validation datasets
* `TCRseq_script_12.R`: Compare abundance of metaclones and published Mtb-reactive TCRs
* `TCRseq_script_13.R`: Compare abundance of metaclones and private Mtb-reactive TCRs
* `TCRseq_script_14.R`: Assess publicity of Mtb-reactive TCRs in day 7 TST

**Notes:**
* `TCRseq_script_4` and `TCRseq_script_5` were adapted from code by Dr Andreas Tiffeau-Mayer. Change `chain` and `mhc_class` in repeated runs of the `HLA dependence` chunk of these scripts to test separately MHC II and MHC I association of cross-donor convergence.
* `TCRseq_script_8` was adapted from code by Prof Benny Chain. This script also produces **Figure S7A**.
* `TCRseq_script_10` identifies the index of selected metaclones. Their cluster number can be looked up in Table S4. The associated adjacency graphs and sequence logos, identifiable by metaclone cluster number, can then be found in the output produced by `Script_5` from the metaclone directory of this repository, and are displayed in **Figure 4D**.

**To make main figures:**
* **Figure 2A**: TCRseq_script_1 &rarr; TCRseq_script_2 &rarr; Plotting_script_Figure2
* **Figure 2B**: TCRseq_script_1 &rarr; TCRseq_script_3 &rarr; Plotting_script_Figure2
* **Figure 2C-E**: TCRseq_script_1 &rarr; TCRseq_script_4 &rarr; Plotting_script_Figure2
* **Figure 3**: TCRseq_script_1 &rarr; TCRseq_script_6 &rarr; TCRseq_script_7 &rarr; Plotting_script_Figure3
* **Figure 5A**: TCRseq_script_1 &rarr; TCRseq_script_11 &rarr; Plotting_script_Figure5
* **Figure 5B**: TCRseq_script_1 &rarr; TCRseq_script_9 &rarr; Plotting_script_Figure5 
* **Figure 5C**: TCRseq_script_1 &rarr; TCRseq_script_6 and TCRseq_script_7 &rarr; TCRseq_script_13 &rarr; Plotting_script_Figure5
* **Figure 5D**: TCRseq_script_1 &rarr; TCRseq_script_3 and TCRseq_script_9 &rarr; TCRseq_script_14 &rarr; Plotting_script_Figure5

**To make supplementary figures:**
* **Figure S3**: TCRseq_script_1 &rarr; TCRseq_script_2 &rarr; Plotting_script_FigureS3
* **Figure S4**: TCRseq_script_1 &rarr; TCRseq_script_3 &rarr; Plotting_script_FigureS4
* **Figure S5**: TCRseq_script_1 &rarr; TCRseq_script_4 and TCRseq_script_5 &rarr; Plotting_script_FigureS5
* **Figure S6**: TCRseq_script_1 &rarr; TCRseq_script 6 &rarr; TCRseq_script_7 &rarr; Plotting_script_FigureS6
* **Figure S7B**: TCRseq_script_1 &rarr; TCRseq_script 6 &rarr; TCRseq_script_8 &rarr; Plotting_script_FigureS7B
* **Figure S8**: TCRseq_script_1 &rarr; TCRseq_script_9 &rarr; Plotting_script_FigureS8 
* **Figure S9**: TCRseq_script_1 &rarr; TCRseq_script_12 &rarr; Plotting_script_FigureS9
* **Figure S10**: TCRseq_script_1 &rarr; TCRseq_script_6 and TCRseq_script_7 &rarr; TCRseq_script_13 &rarr; Plotting_script_FigureS10
* **Figure S11**: TCRseq_script_1 &rarr; TCRseq_script_3 and TCRseq_script_9 &rarr; TCRseq_script_14 &rarr; Plotting_script_FigureS11
