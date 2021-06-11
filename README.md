<!-- <p align="center">
  <img src="LRBinner_logo.png" width="600" title="Final Labelling" alt="Final Labelling">
</p> -->

# LRBinner: Binning Error-Prone Long Reads Using Auto Encoders

![GitHub](https://img.shields.io/github/license/anuradhawick/LRBinner)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/anuradhawick/LRBinner)

# :exclamation:Notice :stop_sign:

* Training on CUDA is now available!

# :exclamation:Notice :stop_sign:

* A new update will be available end of June 2021 with much faster vectorization and GPU support. Stay tuned 😄
* New test examples and results with blogs will also be available soon 🔖 
* We will put a new release as well 👌

## Dependencies
LRBinner is coded purely using C++ (v9) and Python 3.7. To run LRBinner, you will need to install the following python and C++ modules.

### Python dependencies
* numpy 1.16.4 
* scipy 1.3.0 
* seaborn 0.9.0
* h5py 2.9.0
* tabulate 0.8.7
* pytorch 1.4.0

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
* To install the program 
```
pip install .
```
OR add the program path to your $PATH variable.

### Test run data

Extract test data from [here](https://anu365-my.sharepoint.com/:f:/g/personal/u6776114_anu_edu_au/EnV-rUq01pRHl1lH4Y8SaSwBwVVMKNAptbA6YW8RWX6Pqw?e=tDgy9v);

```
LRBinner -r reads.fasta -o lrb_output/ --ae-epochs 200 --resume -mbs 1000 -bit 0 -bs 10 -bc 10
```

### Test run results

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

### Available LRBinner Commands 

Use the `-h` argument to list all the available commands.
```
cd LRBinner
./LRBinner -h
```
### Help

```
usage: LRBinner [-h] --reads-path READS_PATH [--k-size {3,4,5}]
                [--bin-size BIN_SIZE] [--bin-count BIN_COUNT]
                [--ae-epochs AE_EPOCHS] [--ae-dims AE_DIMS]
                [--ae-hidden AE_HIDDEN] [--min-bin-size MIN_BIN_SIZE]
                [--bin-iterations BIN_ITERATIONS] [--threads THREADS]
                [--separate-reads] [--cuda] [--resume] --output <DEST>
                [--version]

LRBinner Help. A tool developed for binning of metagenomics long reads
(PacBio/ONT). Tool utilizes composition and coverage profiles of reads based
on k-mer frequencies to perform dimension reduction via a deep variational
auto-encoder. Dimension reduced reads are then clustered using a novel
distance histogram based clustering algorithm. Minimum RAM requirement is 9GB.

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
  --min-bin-size MIN_BIN_SIZE, -mbs MIN_BIN_SIZE
                        The minimum number of reads a bin should have.
  --bin-iterations BIN_ITERATIONS, -bit BIN_ITERATIONS
                        Number of iterations for cluster search. Use 0 for
                        exhaustive search.
  --threads THREADS, -t THREADS
                        Thread count for computations
  --separate-reads, -sep
                        Flag to separate reads into bins detected. Avaialbe in
                        folder named 'binned'.
  --cuda                Whether to use CUDA if available.
  --resume              Continue from the last step or the binning step (which
                        ever comes first). Can save time needed count k-mers.
  --output <DEST>, -o <DEST>
                        Output directory
  --version, -v         Show version.

```
* Output path is the foldername that you wish the results to be in.
* Specify the number of threads
<!-- 
## Citation

```
TBD
``` -->

## :exclamation:Note :stop_sign:

CODE IS UNDER CLEANING! CHANGES WILL FOLLOW
