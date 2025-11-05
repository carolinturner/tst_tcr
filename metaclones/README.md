# Scripts for metaclone discovery analyses

The scripts provided here were adapted from code written by Dr Andreas Mayer-Tiffeau, and used for metaclone discovery analyses in the TST TCR manuscript, to identify metaclonotypist- and GLIPH2-derived T cell metaclones. The output of these scripts is required as input for some of the _TCRseq_scripts_ and _Plotting_scripts_ listed in the TCRseq directory of this GitHub repository. Relevant output tables are also summarised in Supplementary Files S5-S8 of the manuscript. 

At time of peer-reviewed publication of the manuscript, the processed TCRseq data and associated metadata that are needed to run the scripts can be downloaded from [UCL's Research Data Repository](add here link to DOI 10.5522/04/28049606 once public). HLA imputations for participants are available as Supplementary File S4 of the manuscript. All source data should be saved in a subdirectory called `data`.
 
Metadata = `metadata.csv`
HLA data = Supplementary File S4 saved as `FileS4.csv`

**Overview of analysis scripts:**
* `Metaclones_script_1.R`: Re-format HLA imputation data
* `Metaclones_script_2.ipynb`: Prepare data for metaclone discovery (including down-sampling)
* `Metaclones_script_3.py` and `Snakefile`: Parameter sweep to identify optimal Metaclonotypist parameters for this dataset
* `Metaclones_script_4.ipynb`: Tidy up and explore output from parameter sweep through plots
* `Metaclones_script_5.ipynb`: Metaclonotypist analysis
* `Metaclones_script_6.ipynb`: Gliph2 analysis, followed by HLA association test as implemented in Metaclonotypist
* `Metaclones_script_7.R`: Add regex pattern to Gliph2 clusters
* `Metaclones_script_8.R`: Check metaclone CDR3s against VDJdb Mtb TCRs
* `Metaclones_script_9.ipynb`: Applying different filters to Gliph2 output

**Notes:**
* Down-sampling in `Metaclones_script_2` has been done without setting a reproducibility seed. Downstream outputs may therefore differ slightly from our published analysis.
* `Metaclones_script_3` is run as snakemake workflow, for example with the following terminal command:
	````console
	$ snakemake -s Snakefile -c 8
	````
	Output files are saved in a directory called `results`.
* `Metaclones_script_4` combines all output files from the `results` folder into one data frame.
* Change `MHC class` in repeated runs of `Metaclones_script_5` to test separately for MHC II or MHC I associations in metaclonotypist analysis. 
* Repeated analysis of the same dataset using metaclonotypist (`Metaclones_script_5`) with identical settings may result in slightly different results. This is due to stochasticity of the Leiden clustering step of the metaclonotypist function, implemented through `igraph`, where no reproducibility seed can be supplied.
* Install Gliph2 in a Linux or Mac environment (http://50.255.35.37:8080/tools) before running `Script_6`
* Change `mhc_class` in repeated runs of `Metaclones_script_6` to test separately MHC II or MHC I associations of Gliph2 patterns.
* Change `mhc_class` in repeated runs of `Metaclones_script_8` to process separately Gliph2 outputs of MHC II and MHC I associated gliph patterns.

**To make:**
* **Figure 4C-D**: Metaclones_script_1 &rarr; Metaclones_script_2 &rarr; Metaclones_script_3 and Snakefile &rarr; Metaclones_script_5 &rarr; Plotting_script_Figure4CD
* **Figure 4F**: Metaclones_script_8 &rarr; [deepVenn](https://www.deepvenn.com/)
* **Supplementary Figure 9A**: Metaclones_script_8
* **Table 3**: Metaclones_Table_script