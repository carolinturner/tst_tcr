# Scripts for TCRseq analyses

The scripts provided here were used for TCRseq analyses in the [TST TCR publication](link here once available). _TCRseq_scripts_ generate the data required for _Plotting_scripts_ to make Figures 2-5 and Supplementary Figures S2-S11

**Input data:**
* Processed TCRseq data and metadata can be downloaded from [UCL's Research Data Repository](add here link to DOI 10.5522/04/28049606 once public).
* Antigen-reactive CDR3 sequences collated from published sources are available as Table S2 of the [manuscript](link here once available).
* HLA imputations for participants are available as Table S3.
* Metaclone summaries (including HLA association and regular expression patterns) are available as Tables S4-S9. These metaclones were discovered using the code availabe in the metaclone directory of this GitHub repository.
* All source data should be saved in a subdirectory called `data`.
* Supplementary tables should be saved as individual .csv files, e.g. `TableS2.csv`

**Overview of scripts:**
* `TCRseq_script_1.R`: Pre-processing of data
* `TCRseq_script_2.R`: Calculate diversity indices
* `TCRseq_script_3.R`: Calculate abundance of published antigen-reactive TCRs
* `TCRseq_script_4.ipynb`: Coincidence analysis (down-sampled repertoires)
* `TCRseq_script_5.ipynb`: Coincidence analysis (full repertoires)
* `TCRseq_script_6.R`: Define private antigen-reactive TCRs
* `TCRseq_script_7.R`: Calculate abundance of private antigen-reactive TCRs
* `TCRseq_script_8.R`: Expansion of antigen-reactive TCRs in paired TST samples from day 2 and day 7

**Notes:**
* `TCRseq_script_4` and `TCRseq_script_5` were adapted from code by Dr Andreas Tiffeau-Mayer. Change `chain` and `mhc_class` in repeated runs of the `HLA dependence` chunk of these scripts to test separately MHC II and MHC I association of cross-donor convergence.
* `TCRseq_script_8` was adapted from code by Prof Benny Chain. This script also produces Figure S6A.

**To make main figures:**
* **Figure 2A**: TCRseq_script_1 &rarr; TCRseq_script_2 &rarr; Plotting_script_Figure2A
* **Figure 2B**: TCRseq_script_1 &rarr; TCRseq_script_3 &rarr; Plotting_script_Figure2B
* **Figure 2C-E**: TCRseq_script_1 &rarr; TCRseq_script_4 &rarr; Plotting_script_Figure2C-E
* **Figure 3**: TCRseq_script_1 &rarr; TCRseq_script_6 and TCRseq_script_7 &rarr; Plotting_script_Figure3

**To make supplementary figures:**
* **Figure S2**: TCRseq_script_1 &rarr; TCRseq_script_2 &rarr; Plotting_script_FigureS2
* **Figure S3**: TCRseq_script_1 &rarr; TCRseq_script_3 &rarr; Plotting_script_FigureS3
* **Figure S4**: TCRseq_script_1 &rarr; TCRseq_script_4 and TCRseq_script_5 &rarr; Plotting_script_FigureS4
* **Figure S5**: TCRseq_script_1 &rarr; TCRseq_script 6 and TCRseq_script_7 &rarr; Plotting_script_FigureS5
* **Figure S6**: TCRseq_script_1 &rarr; TCRseq_script 6 and TCRseq_script_8 &rarr; Plotting_script_FigureS6

