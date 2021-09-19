<!-- <p align="center">
  <img src="LRBinner_logo.png" width="600" title="Final Labelling" alt="Final Labelling">
</p> -->

# LRBinner: Binning Error-Prone Long Reads Using Auto Encoders

![GitHub](https://img.shields.io/github/license/anuradhawick/LRBinner)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/anuradhawick/LRBinner)

## Updates

<!-- * Reference based binning is now available (check Usage section) -->
<!-- * Training on CUDA is now available! -->
* Contigs from long reads assemblies can now be binned! 
* Taxonomic binning of long reads paused till next update

## Dependencies
LRBinner is coded purely using C++ (v9) and Python 3.7. To run LRBinner, you will need to install the following python and C++ modules.

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

Extract test data from [here](https://anu365-my.sharepoint.com/:f:/g/personal/u6776114_anu_edu_au/EnV-rUq01pRHl1lH4Y8SaSwBwVVMKNAptbA6YW8RWX6Pqw?e=tDgy9v);

```
LRBinner -r reads.fasta -o lrb_output/ --ae-epochs 200 --resume -mbs 1000 -bit 0 -bs 10 -bc 10
```

## Test run results

Note that the results could vary due to slight sampling differences. Evaluations can be done using the `eval.py` script.

```
_                         Bin-0_(2)  Bin-1_(3)  Bin-2_(1)  Bin-3_(0)  Bin-4_(4)
Campylobacter_jejuni      1          0          0          63         5726
Acetobacter_pasteurianus  592        6138       58         42         0
Yersinia_pestis_Angola    30946      255        25         370        181
Lactobacillus_paracasei   157        0          6962       100        0
Chlamydia_psittaci        0          0          0          5512       0

Precision            96.77
Recall               96.77
F1-Score             96.77
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

> More updates to come!

> Get in touch [anuradhawick.com](https://anuradhawick.com)
