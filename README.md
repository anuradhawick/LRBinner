<!-- <p align="center">
  <img src="MetaBCC-LR_logo.png" width="600" title="Final Labelling" alt="Final Labelling">
</p> -->

# MetaBCC-LR-2: Binning Error-Prone Long Reads Using Auto Encoders

![GitHub](https://img.shields.io/github/license/anuradhawick/MetaBCC-LR-2)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/anuradhawick/MetaBCC-LR-2)

## Dependencies
MetaBCC-LR-2 is coded purely using C++ (v9) and Python 3.7. To run MetaBCC-LR-2, you will need to install the following python and C++ modules.

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

### Third party programs
* DSK: https://github.com/GATB/dsk
    * Add DSK binaries to your PATH variable

## Downloading MetaBCC-LR
To download MetaBCC-LR-2, you have to clone the MetaBCC-LR-2 repository to your machine.

```
git clone https://github.com/anuradhawick/MetaBCC-LR-2.git
```

## Compiling the source code
* Build the binaries
```
cd MetaBCC-LR-2
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

## Running the MetaBCC-LR
In order to run MetaBCC-LR you are required to provide the reads in FASTQ or FASTA format.

```
cd MetaBCC-LR-2
./MetaBCC-LR-2 -h


```
* Output path is the foldername that you wish the results to be in.
* Specify the number of threads
* The program requires a minimum of 5GB to run. This is because we have optimized the coverage histogram generation process to accommodate all 15mers in RAM for faster lookup of counts.
<!-- 
## Citation

```
TBD
``` -->

## Notes

CODE IS UNDER CLEANING! CHANGES WILL FOLLOW
