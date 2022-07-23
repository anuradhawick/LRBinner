<!-- <p align="center">
  <img src="LRBinner_logo.png" width="600" title="Final Labelling" alt="Final Labelling">
</p> -->

# LRBinner: Binning Error-Prone Long Reads Using Auto Encoders

![GitHub](https://img.shields.io/github/license/anuradhawick/LRBinner)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/anuradhawick/LRBinner)

## Use Docker!

`Dockerfile` is now available. (If you're familiar deploying a docker file, No Image pushed yet)

## Dependencies
LRBinner is coded purely using C++ (v9) and Python 3.7. To run LRBinner, you will need to install the following python and C++ modules.

A possible conda environment to work (credits [Calum Walsh
](https://github.com/cazzlewazzle89))

```sh
conda create -n lrbinner -y python=3.7 numpy scipy seaborn h5py tabulate pytorch hdbscan gcc openmp tqdm biopython fraggenescan hmmer
```

### Python dependencies
Essential libraries

* numpy 1.16.4 
* scipy 1.3.0 
* seaborn 0.9.0
* h5py 2.9.0
* tabulate 0.8.7
* pytorch 1.4.0

Essential for contig binning
<!-- * umap-learn -->
* fraggenescan 1.31
* hmmer 3.3.2
* HDBSCAN


### C++ requirements
* GCC version 9.1.0
* OpenMP 4.5 for multi processing

## Downloading LRBinner
To download LRBinner, you have to clone the LRBinner repository to your machine.

```
git clone https://github.com/anuradhawick/LRBinner.git
```

## Compiling the source code
* Build the binaries
```
cd LRBinner
python setup.py build
```
OR
```
sh build.sh
```    
<!-- * To install the program 
```
pip install .
```
OR add the program path to your $PATH variable. -->

## Test run data 
>For long reads binning coverage of this dataset may be too low for assembly)

Extract Sim-8 data from [here](https://anu365-my.sharepoint.com/:u:/g/personal/u6776114_anu_edu_au/EZfR_olp6fxOpoKEexcm1ZABk1fTQTnW0-3ja772k22WbA?e=aoXr7N);

```
python LRBinner reads -r reads.fasta -bc 10 -bs 32 -o lrb --resume --cuda -mbs 5000 --ae-dims 4 --ae-epochs 200 -bit 0 -t 32
```

## Test run results

Note that the results could vary due to slight sampling differences. Evaluations can be done using the `eval.py` script.

```
_                                                           Bin-0_(4)  Bin-1_(1)  Bin-2_(3)  Bin-3_(0)  Bin-4_(5)  Bin-5_(2)  Bin-6_(6)  Bin-7_(7)
CP002618.1_Lactobacillus_paracasei_strain_BD-II             744        973        92         70169      215        0          0          0
NC_011658.1_Bacillus_cereus_AH187                           110        11         277        10         30129      0          1          439
CP002807.1_Chlamydia_psittaci_08DC60_chromosome             0          0          9          0          1          0          0          11014
NC_012883.1_Thermococcus_sibiricus_MM_739                   74         0          32441      1          28         2          0          8
NC_011415.1_Escherichia_coli_SE11,_complete_sequence        70474      432        47         220        616        3          0          43
NC_013929.1_Streptomyces_scabiei_87.22_complete_genome      343        0          0          0          0          118979     8          0
FQ312002.1_Haemophilus_parainfluenzae_T3T1_complete_genome  394        83623      32         191        1596       0          0          46
AP011170.1_Acetobacter_pasteurianus_IFO_3283-12_DNA         991        40         15         83         14         0          7382       13

Precision	     98.12
Recall    	     98.12
F1-Score  	     98.12
Bins      	         8
```
## Usage
<!-- ### Constraints file format

Provide a TSV file containing read index and the taxonomic label. See the following example. An example usecase would be to sample few reads (10000-50000) and annotate them. Then prepare a TSV file as below and use in binning.

Must-link - reads from same taxa
Cannot-link - reads from different taxa

```
107885  Escherichia_coli
188789  Haemophilus_parainfluenzae
181635  Haemophilus_parainfluenzae
388619  Streptomyces_scabiei
27303   Bacillus_cereus
122292  Escherichia_coli
20565   Bacillus_cereus
99539   Escherichia_coli
300860  Streptomyces_scabiei
294442  Streptomyces_scabiei
```
Use flag `--constraints` or `-c` to provide this file to LRBinner. -->

### Parameters

Our manuscript presents results with `-k 3` i.e. using 3-mers. Use `-k 4` for tetramer based binning.
Internal parameters are not yet set for `-k 5` choice. We are working on that. :smile:

### Available LRBinner Commands 

Use the `-h` argument to list all the available commands.
```
cd LRBinner
./LRBinner -h
```
### Help

```
usage: LRBinner [-h] [--version] {reads,contigs} ...

LRBinner Help. A tool developed for binning of metagenomics long reads
(PacBio/ONT) and long read assemblies. Tool utilizes composition and coverage
profiles of reads based on k-mer frequencies to perform dimension reduction
via a deep variational auto-encoder. Dimension reduced reads are then
clustered. Minimum RAM requirement is 9GB (4GB GPU if cuda used).

optional arguments:
  -h, --help       show this help message and exit
  --version, -v    Show version.

LRBinner running Mode:
  {reads,contigs}
    reads          for binning reads
    contigs        for binning contigs
```
### Reads binning help

```
usage: LRBinner reads [-h] --reads-path READS_PATH [--k-size {3,4,5}]
                      [--bin-size BIN_SIZE] [--bin-count BIN_COUNT]
                      [--ae-epochs AE_EPOCHS] [--ae-dims AE_DIMS]
                      [--ae-hidden AE_HIDDEN] [--threads THREADS] [--separate]
                      [--cuda] [--resume] --output <DEST> [--version]
                      [--min-bin-size MIN_BIN_SIZE]
                      [--bin-iterations BIN_ITERATIONS]

optional arguments:
  -h, --help            show this help message and exit
  --reads-path READS_PATH, -r READS_PATH
                        Reads path for binning
  --k-size {3,4,5}, -k {3,4,5}
                        k value for k-mer frequency vector. Choose between 3
                        and 5.
  --bin-size BIN_SIZE, -bs BIN_SIZE
                        Bin size for the coverage histogram.
  --bin-count BIN_COUNT, -bc BIN_COUNT
                        Number of bins for the coverage histogram.
  --ae-epochs AE_EPOCHS
                        Epochs for the auto_encoder.
  --ae-dims AE_DIMS     Size of the latent dimension.
  --ae-hidden AE_HIDDEN
                        Hidden layer sizes eg: 128,128
  --threads THREADS, -t THREADS
                        Thread count for computations
  --separate, -sep      Flag to separate reads/contigs into bins detected.
                        Avaialbe in folder named 'binned'.
  --cuda                Whether to use CUDA if available.
  --resume              Continue from the last step or the binning step (which
                        ever comes first). Can save time needed count k-mers.
  --output <DEST>, -o <DEST>
                        Output directory
  --version, -v         Show version.
  --min-bin-size MIN_BIN_SIZE, -mbs MIN_BIN_SIZE
                        The minimum number of reads a bin should have.
  --bin-iterations BIN_ITERATIONS, -bit BIN_ITERATIONS
                        Number of iterations for cluster search. Use 0 for
                        exhaustive search.
```
### Contigs binning help

```
usage: LRBinner contigs [-h] --reads-path READS_PATH [--k-size {3,4,5}]
                        [--bin-size BIN_SIZE] [--bin-count BIN_COUNT]
                        [--ae-epochs AE_EPOCHS] [--ae-dims AE_DIMS]
                        [--ae-hidden AE_HIDDEN] [--threads THREADS]
                        [--separate] [--cuda] [--resume] --output <DEST>
                        [--version] --contigs CONTIGS

optional arguments:
  -h, --help            show this help message and exit
  --reads-path READS_PATH, -r READS_PATH
                        Reads path for binning
  --k-size {3,4,5}, -k {3,4,5}
                        k value for k-mer frequency vector. Choose between 3
                        and 5.
  --bin-size BIN_SIZE, -bs BIN_SIZE
                        Bin size for the coverage histogram.
  --bin-count BIN_COUNT, -bc BIN_COUNT
                        Number of bins for the coverage histogram.
  --ae-epochs AE_EPOCHS
                        Epochs for the auto_encoder.
  --ae-dims AE_DIMS     Size of the latent dimension.
  --ae-hidden AE_HIDDEN
                        Hidden layer sizes eg: 128,128
  --threads THREADS, -t THREADS
                        Thread count for computations
  --separate, -sep      Flag to separate reads/contigs into bins detected.
                        Avaialbe in folder named 'binned'.
  --cuda                Whether to use CUDA if available.
  --resume              Continue from the last step or the binning step (which
                        ever comes first). Can save time needed count k-mers.
  --output <DEST>, -o <DEST>
                        Output directory
  --version, -v         Show version.
  --contigs CONTIGS, -c CONTIGS
                        Contigs path
```

## Citation

If you use LRBinner please cite using the the following bibtex files

```bibtex
@InProceedings{wickramarachchi_et_al:LIPIcs.WABI.2021.11,
  author =	{Wickramarachchi, Anuradha and Lin, Yu},
  title =	{{LRBinner: Binning Long Reads in Metagenomics Datasets}},
  booktitle =	{21st International Workshop on Algorithms in Bioinformatics (WABI 2021)},
  pages =	{11:1--11:18},
  series =	{Leibniz International Proceedings in Informatics (LIPIcs)},
  ISBN =	{978-3-95977-200-6},
  ISSN =	{1868-8969},
  year =	{2021},
  volume =	{201},
  editor =	{Carbone, Alessandra and El-Kebir, Mohammed},
  publisher =	{Schloss Dagstuhl -- Leibniz-Zentrum f{\"u}r Informatik},
  address =	{Dagstuhl, Germany},
  URL =		{https://drops.dagstuhl.de/opus/volltexte/2021/14364},
  URN =		{urn:nbn:de:0030-drops-143644},
  doi =		{10.4230/LIPIcs.WABI.2021.11},
  annote =	{Keywords: Metagenomics binning, long reads, machine learning, clustering}
}
```

```bibtex
﻿@Article{Wickramarachchi2022,
  author={Wickramarachchi, Anuradha and Lin, Yu},
  title={Binning long reads in metagenomics datasets using composition and coverage information},
  journal={Algorithms for Molecular Biology},
  year={2022},
  month={Jul},
  day={11},
  volume={17},
  number={1},
  pages={14},
  issn={1748-7188},
  doi={10.1186/s13015-022-00221-z},
  url={https://doi.org/10.1186/s13015-022-00221-z}
}

```

> More updates to come!

> Get in touch [anuradhawick.com](https://anuradhawick.com)
