import os
import sys
from Bio import SeqIO
import pickle
import argparse

parser = argparse.ArgumentParser(description='Separate reads in to bins.')

parser.add_argument('--reads', '-r', help="Path to reads used for binning", type=str, required=True)
parser.add_argument('--bins', '-b', help="Path of binning_result.pkl file from LRBinner", type=str, required=True)
parser.add_argument('--outpath', '-o', help="Output directory (will be created or content inside will be removed before execution)" ,type=str, required=True)

args = parser.parse_args()

reads = args.reads
bins = args.bins
outpath = args.outpath

fmt = "fasta" if reads.split('.')[-1].lower() in ["fasta", "fna", "fa"] else "fastq"
bins = pickle.load(open(bins, "rb"))
bin_files = {}
read_bin = {}

if not os.path.isdir(f"{outpath}/"):
    os.makedirs(f"{outpath}/")

for k, v in bins.items():
    if not os.path.isdir(f"{outpath}/bin-{k}/"):
        os.makedirs(f"{outpath}/bin-{k}/")

    bin_files[k] = open(f"{outpath}/bin-{k}/reads.fasta", "w+")

    for r in v:
        read_bin[r] = k

binout = open(f"{outpath}/bins.txt", "w+")
lenout = open(f"{outpath}/lengths.txt", "w+")

for r, record in enumerate(SeqIO.parse(reads, fmt)):
    binout.write(f"{read_bin[r]}\n")
    lenout.write(f"{len(record.seq)}\n")
    bin_files[read_bin[r]].write(f">read-{r}\n")
    bin_files[read_bin[r]].write(f"{record.seq}\n")

binout.close()
lenout.close()

for k, f in bin_files.items():
    f.close()