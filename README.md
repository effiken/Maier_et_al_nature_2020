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
   - [tglkmeans](https://github.com/tanaylab/tglkmeans)
   - [scDissector](https://github.com/effiken/scDissector)
3. Downloaded and unzipped version of this repository  on a local path.

Running the scripts in R
Assuming Maier_et_al_nature_2020 is the local path of the repository, the following script will run the clustering distributedly on LSF:

source("Maier_et_al_nature_2020/scripts/clustering/run_clustering_mouse.r")

Note: Each run of the clustering might produce slightly different results due to different random seeds.
