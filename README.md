# Pseudotyped lentivirus neutralization assays

This repo includes all data from the pseudotyped lentivirus neutralization assays for [Ellis, et al. 2021](). 

## Repo organization

The [MouseNeuts.ipynb](https://github.com/jbloomlab/RBD_nanoparticle_vaccine/blob/organize_repo/MouseNeuts.ipynb) notebook contains analyses using [neutcurve](https://jbloomlab.github.io/neutcurve/) to calculate IC50 values for serum samples from mice immunized with RBD vaccine constructs as described in the manuscriptlinked above. 

The [all_neut_results.csv](https://github.com/jbloomlab/RBD_nanoparticle_vaccine/blob/organize_repo/all_neut_results.csv) file contains IC50s and additioonal metadata for all samples run for this project, including samples that were run more than once.

The [mouse_neuts.csv](https://github.com/jbloomlab/RBD_nanoparticle_vaccine/blob/organize_repo/mouse_neuts.csv) file contains IC50s and additional metadata for the mouse serum samples. This file only includes data for the final run of the samples that were run more than once. These data are presented in Figure 4 of the manuscript.

The [mouse_plus_ctrls_neuts.csv](https://github.com/jbloomlab/RBD_nanoparticle_vaccine/blob/organize_repo/mouse_plus_ctrls_neuts.csv) file contains only the final data for the mouse neuts and all IC50s for the control samples. 

### Description of directories

* [MouseNeuts_files](https://github.com/jbloomlab/RBD_nanoparticle_vaccine/tree/organize_repo/MouseNeuts_files): Output files from markdown conversion of `MouseNeuts.ipynb` notebook.

* [RLU_data](https://github.com/jbloomlab/RBD_nanoparticle_vaccine/tree/organize_repo/RLU_data): Excel files from the Bloom lab plate reader with raw luciferase data for each assay.

* [fract_infect](https://github.com/jbloomlab/RBD_nanoparticle_vaccine/tree/organize_repo/fract_infect): CSV files with calculated fraction infectivity values for each assay. The `excel_to_fractinfect.py` script was used to convert the RLU data and sample maps into the fraction infectivity csv files.

* [sample_maps](https://github.com/jbloomlab/RBD_nanoparticle_vaccine/tree/organize_repo/sample_maps): Sample info needed to use `excel_to_fractinfect.py` to convert RLU data to fraction infectivities.

* [PlateLayouts](https://github.com/jbloomlab/RBD_nanoparticle_vaccine/tree/organize_repo/PlateLayouts): Plate layouts for neutralization assays.

 
