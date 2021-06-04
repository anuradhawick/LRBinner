import pickle
import os
import logging
import sys

logger = logging.getLogger('LRBinner')


def run_kmers(reads_path, output, k_size, threads):
    if not os.path.isdir(f"{output}/profiles"):
        os.makedirs(f"{output}/profiles")

    cmd = f""""{os.path.dirname(__file__)}/bin/count-kmers" "{reads_path}" "{output}/profiles/com_profs" {k_size} {threads}"""
    logger.debug("CMD::" + cmd)
    o = os.system(cmd)
    check_proc(o, "Counting Trimers")


def run_15mer_counts(reads_path, output, threads):
    if not os.path.isdir(f"{output}/profiles"):
        os.makedirs(f"{output}/profiles")

    cmd = f""""{os.path.dirname(__file__)}/bin/count-15mers" "{reads_path}" "{output}/profiles/15mers-counts" {threads}"""
    logger.debug("CMD::" + cmd)
    o = os.system(cmd)
    check_proc(o, "Counting 15-mers")


def run_15mer_vecs(reads_path, output, bin_size, bins, threads):
    if not os.path.isdir(f"{output}/profiles"):
        os.makedirs(f"{output}/profiles")

    cmd = f""""{os.path.dirname(__file__)}/bin/search-15mers" "{output}/profiles/15mers-counts" "{reads_path}" "{output}/profiles/cov_profs" {bin_size} {bins} {threads}"""
    logger.debug("CMD::" + cmd)
    o = os.system(cmd)
    check_proc(o, "Counting 15-mer profiles")


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
        if name != "":
            logger.error(f"Error in step: {name}")
        logger.error("Failed due to an error. Please check the log. Good Bye!")
        sys.exit(ret)
