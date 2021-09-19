import pickle
import os
import logging
import sys
from collections import defaultdict
from Bio import SeqIO
import numpy as np

from metacoag_utils import marker_gene_utils
from mbcclr_utils import ae_utils
from mbcclr_utils import cluster_utils

logger = logging.getLogger('LRBinner')


class Checkpointer():
    def __init__(self, checkpoint_path, _load_to_resume=False):
        self.cpath = checkpoint_path
        if _load_to_resume and os.path.isfile(self.cpath):
            self.completed = pickle.load(open(self.cpath, "rb"))
        else:
            self.completed = {}


    def should_run_step(self, stage, params):
        if stage not in self.completed:
            return True
        return self.completed[stage] != params


    def log(self, stage, params):
        self.completed[stage] = params
        ps, pc = [int(x) for x in stage.split("_")][:2]

        # remove all child stages
        for s in list(self.completed.keys()):
            p, c = [int(x) for x in s.split("_")][:2]
            # is stage is seen after current checkpoint
            if p > ps:
                del self.completed[s]

        self._save()


    def _save(self):
        pickle.dump(self.completed, open(self.cpath, "wb+"))


    def __str__(self):
        return str(self.completed)


def split_contigs(contigs, output):
    contig_groups = defaultdict(list)
    fragment_parent = {}

    with open(f"{output}/fragments/contigs.fasta", "w+") as scf:
        i = 0
        for n, record in enumerate(SeqIO.parse(contigs, "fasta")):
            if len(record.seq) >= 5000:
                sub_contigs = [record.seq[x:x+2500] for x in range(0, len(record.seq), 2500)]
                sub_contigs.append(record.seq[-2500:])
            else:
                sub_contigs = [record.seq]

            for sc in sub_contigs:
                rid = f">{n}_{i}"
                rec = str(sc)
                scf.write(f"{rid}\n{rec}\n")

                contig_groups[str(record.id)].append(i)
                fragment_parent[i] = str(record.id)
                i += 1                

    return contig_groups, fragment_parent


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


def run_15mer_vecs(reads_path, output, bin_size, bin_count, threads):
    if not os.path.isdir(f"{output}/profiles"):
        os.makedirs(f"{output}/profiles")

    cmd = f""""{os.path.dirname(__file__)}/bin/search-15mers" "{output}/profiles/15mers-counts" "{reads_path}" "{output}/profiles/cov_profs" {bin_size} {bin_count} {threads}"""
    logger.debug("CMD::" + cmd)
    o = os.system(cmd)
    check_proc(o, "Counting 15-mer profiles")


def check_proc(ret, name=""):
    if ret != 0:
        if name != "":
            logger.error(f"Error in step: {name}")
        logger.error("Failed due to an error. Please check the log. Good Bye!")
        sys.exit(ret)
