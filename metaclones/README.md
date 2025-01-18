# Scripts for metaclone discovery analyses

The Jupyter notebooks provided here were adapted from code written by Andreas Mayer-Tiffeau, and used for metaclone discovery analyses in the [TST TCR publication](link here once available), to identify metaclonotypist- and gliph2-derived T cell metaclones. The output of these scripts is required as input for some of the _TCRseq_scripts_ and _Plotting_scripts_ listed in the TCRseq directory of this GitHub repository. Most output tables are also summarised in Supplementary Tables S3-S8 of the [manuscript](link here once available).

To run the code, processed TCRseq data and associated metadata can be downloaded from [UCL's Research Data Repository](add here link to DOI 10.5522/04/28049606 once public). HLA imputations for participants are available as Table S2 of the [manuscript](link here once available). All source data should be saved in a subdirectory called `data`.
 
Metadata = `metadata.csv`
HLA data = Table S2 saved `TableS2.csv`

**Overview of scripts:**
* `Script_1.R`: Re-format HLA imputation data
* `Script_2.ipynb`: Prepare data for metaclone discovery (including down-sampling)
* `Script_3.ipynb`: Metaclonotypist analysis
* `Script_4.ipynb`: Gliph2 analysis
* `Script_5.R`: Tidy up Gliph2 output (including conversion to regex patterns)

_Notes:_
* Down-sampling in `Script_2` has been done without setting a reproducibility seed. Downstream outputs may therefore differ slightly from our published analysis.
* Change `chain` and `MHC class` in repeated runs of `Script_3`to perform metaclonotypist analysis separately on alpha or beta chain repertoires, testing separately for MHC II or MHC I associations. 
* Repeated analysis of the same dataset using metaclonotypist (`Script_3`) with identical settings may result in slightly different results. This is due to stochasticity of the Leiden clustering step of the metaclonotypist function, implemented through `igraph`, where no reproducibility seed can be supplied.
* Change `MHC class` in repeated runs of `Script_4` to test separately MHC II or MHC I associations of Gliph2 patterns (applied to beta chain repertoires only).
* Change `mhc_class` in repeated runs of `Script_5` to process separately Gliph2 outputs of MHC II and MHC I associated gliph patterns.