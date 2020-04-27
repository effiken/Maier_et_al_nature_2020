
## MAIER_ET_AL_NATURE_2020

This repository contains scripts used for analysis and figure generation in Maier, B., Leader, A.M., Chen, S.T. et al. A conserved dendritic-cell regulatory program limits antitumour immunity. Nature 580, 257–262 (2020). https://doi.org/10.1038/s41586-020-2134-y

NOTE: Clustering and analysis of associated human lung cancer dataset will be reported in partner publication Leader, et al. (in preparation as of 4/27/20). 

As of 4/27/20, this repository is being adapted so that users can download and immediately reproduce all figures locally. Please bear with us, and contact andrew.leader@icahn.mssm.edu if you have any questions.

All mice sequencing data are publicly available (GEO accession code GSE131957). All human sequencing data is available on NCBI with BioProject ID PRJNA609924.

## Making the figures
### Requirements

Tested on Windows 10

1. R
2. R packages: 
	-gplots
	-Matrix.utils
	-mixtools
	-seriation
	-sp
	-scales
	-[scDissector](https://github.com/effiken/scDissector)

3. Downloaded and unzipped version of this repository  on a local path.

### Running the scripts in R

Assuming Maier_et_al_nature_2020 is the local path of the repository we need to load the script files:

source("plot_figures/figures_main.R")`

### Output

Figure will be generated in:
  - plot_figures/

## Clustering

### Requirements

Tested on linux LSF HPC. Due to lack of support of some of the depdendencies, the script cannot run on macOS or Windows.

1. R
2. R packages:
   - Matrix
   - Matrix.utils
   - gplots
   - seriation
   - [tglkmeans](https://bitbucket.org/tanaylab/tglkmeans)
   - [scDissector](https://github.com/effiken/scDissector)
3. Downloaded and unzipped version of this repository  on a local path.

### Running the scripts in R

Assuming Maier_et_al_nature_2020 is the local path of the repository, the following script will run the clustering distributedly on LSF:

`source("Maier_et_al_nature_2020/clustering/run_clustering_[].r")`

Note: Each run of the clustering might produce slightly different results due to different random seeds.
