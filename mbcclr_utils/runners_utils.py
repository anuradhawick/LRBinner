import pickle
import os
from Bio import SeqIO
from tqdm import tqdm
import logging
import sys

from . import scan_dsk

logger = logging.getLogger('LRBinner')

def run_kmers(reads, output, k_size, threads):
    if not os.path.isdir(os.path.dirname(output)):
        os.makedirs(os.path.dirname(output))

    cmd = f""""{os.path.dirname(__file__)}/bin/count-kmers" "{reads}" "{output}" {k_size} {threads} """
    logger.debug(cmd)
    o = os.system(cmd)
    check_proc(o, "Counting Trimers")

def run_15mers(dsk_file, reads, output, bin_size, bin_count, threads):
    if not os.path.isdir(os.path.dirname(output)):
        os.makedirs(os.path.dirname(output))

    cmd = f""""{os.path.dirname(__file__)}/bin/search-15mers" "{dsk_file}" "{reads}" "{output}" {bin_size} {bin_count} {threads}"""
    logger.debug(cmd)
    o = os.system(cmd)
    check_proc(o, "Counting 15-mer profiles")

def run_dsk(reads, output, min_abundance, max_memory, threads):
    if not os.path.isdir(f"{output}"):
        os.makedirs(f"{output}")
    logger.debug("Running DSK")
    cmdDSK = f"""dsk -verbose 0 -file "{reads}" -kmer-size 15 -abundance-min {min_abundance} -out-dir "{output}/DSK" -max-memory {max_memory} -nb-cores {threads}"""
    logger.debug(cmdDSK)
    o = os.system(cmdDSK)
    check_proc(o, "Running DSK")
    o = os.system(f"""mv "{output}"/DSK/*.h5 "{output}"/DSK/output.h5 """)
    check_proc(o, "Rename DSK")
    scan_dsk.scan_dsk(f"{output}/DSK/output.h5", threads, f"{output}/DSK/")

def checkpoint(new_checkpoints, path):
    pickle.dump(new_checkpoints, open(path, "wb+"))

def load_checkpoints(path):
    if not os.path.isfile(path):
        data = {}
        data['completed'] = set()
        return data
    else:
        return pickle.load(open(path, "rb"))

def check_proc(ret, name=""):
    if ret != 0:
        if name!= "": logger.error(f"Error in step: {name}")
        logger.error("Failed due to an error. Please check the log. Good Bye!")
        sys.exit(ret)
