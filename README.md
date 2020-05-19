NOTE: Clustering and analysis of associated human lung cancer dataset will be reported in partner publication Leader, et al. (submitted as of 5/19/20). 

All mice sequencing data are publicly available (GEO accession code GSE131957). All human sequencing data is available on NCBI with BioProject ID PRJNA609924.

## Downloading the data

Data can either be downloaded automatically by running the script to reproduce the figures (see below).
Alternatively, .rd files can be downloaded using the following dropbox links:

mouse scRNAseq data: https://www.dropbox.com/s/3ffq3z75af37rgr/mouse_dc.rd?dl=1

human DC scRNA & CITEseq data: https://www.dropbox.com/s/e0cvefqzoy4hz0r/human_dc.rd?dl=1

These links will download R data structures, the components of which contain the count matrices and cell metadata.

For the human DC dataset, the data structure is called "human_dc" and has the following components:
1. filtered_umitab is the counts table
2. cell_to_tissue, cell_to_sample, cell_to_annot, cell_to_patient are arrays with metadata identifying the tissue, sample, and patient of origin as well as our annotation (cDC1, cDC2, mregDC)
3. filtered_ds is a downsampled matrix, where we have randomly drawn 2000 UMI from each cell containing at least that many. We use this sampling of UMI as a normalization for visualization that doesn't cause artifacts.
4. adt_matrix_by_sample is a list of CITEseq ADT matrices for the samples that had CITEseq. We recommend within-patient normalization for analyses.


## Making the figures
### Requirements

Tested on Windows 10

1. R
2. R packages: 
	- gplots
	- Matrix.utils
	- mixtools
	- seriation
	- sp
	- scales
	- [scDissector](https://github.com/effiken/scDissector)

3. Downloaded and unzipped version of this repository  on a local path.

### Running the scripts in R

Assuming Maier_et_al_nature_2020 is the local path of the repository we need to load the script files:

source("figures_main.R")`

### Output

The above referenced dropbox links will download automatically to a new /data/ directory.

Figure will be generated in a new directory:
  - figures_out/

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

source("clustering/run_clustering_mouse.r")

Note: Each run of the clustering might produce slightly different results due to different random seeds.
