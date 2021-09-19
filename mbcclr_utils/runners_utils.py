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


# stage in form step_substep, dict has step params as values array
def mark_checkpoint(new_checkpoints, path, stage=None):
    if stage is None:
        pickle.dump(new_checkpoints, open(path, "wb+"))
    else:
        ps, pc = [int(x) for x in stage.split("_")]
        saveable = {}
        # add parallel tasks to saveable
        # skip child tasks
        for s in list(new_checkpoints.keys()):
            p, c = [int(x) for x in s.splti("_")]
            if ps >= p:
                saveable[s] = new_checkpoints[s]

        pickle.dump(saveable, open(path, "wb+"))


def should_run_step(stage_params, old_checkpoint, stage):
    if stage not in old_checkpoint:
        return True
    return stage_params != old_checkpoint[stage]

def load_checkpoints(path):
    if not os.path.isfile(path):
        data = {}
        return data
    else:
        return pickle.load(open(path, "rb"))


def check_proc(ret, name=""):
    if ret != 0:
        if name != "":
            logger.error(f"Error in step: {name}")
        logger.error("Failed due to an error. Please check the log. Good Bye!")
        sys.exit(ret)


def run_contig_binning(args):
    # commong arguments for reads binning and contig binning
    reads_path = args.reads_path
    threads = args.threads
    bin_size = args.bin_size
    bin_count = args.bin_count
    k_size = args.k_size
    epochs = args.ae_epochs
    dims = args.ae_dims
    hidden = list(map(int, args.ae_hidden.split(",")))
    binreads = args.separate_reads
    cuda = args.cuda
    resume = args.resume
    min_cluster_size = max(args.min_bin_size, 1)
    contigs = args.contigs
    output = args.output

    checkpoints_path = f"{output}/checkpoints"

    if not resume:
        checkpoint = {}
        mark_checkpoint(checkpoint, checkpoints_path)
    else:
        logger.info("Resuming the program from previous checkpoints")
        checkpoint = load_checkpoints(checkpoints_path)
        logger.debug(str(checkpoint))

    stage = "1_1"
    stage_params = ['contigs_binning']

    if should_run_step(stage_params, checkpoint, stage):
        checkpoint[stage] = stage_params
        mark_checkpoint(checkpoint, checkpoints_path)

    # compute contig lengths
    stage = "2_1"
    stage_params = [contigs]

    if should_run_step(stage_params, checkpoint, stage):
        logger.info("Computing contig lengths")
        contig_length = {}

        for record in SeqIO.parse(contigs, "fasta"):
            contig_length[record.id] = len(record.seq)
        
        pickle.dump(contig_length, open(f"{output}/profiles/contig_lengths.pkl", "wb+"))

        checkpoint[stage] = stage_params
        mark_checkpoint(checkpoint, checkpoints_path)
    else:
        logger.info("Loading contig lengths")
        contig_length = pickle.load(open(f"{output}/profiles/contig_lengths.pkl", "rb"))
    
    # splitting contigs
    stage = "2_2"
    stage_params = [contigs]

    if should_run_step(stage_params, checkpoint, stage):
        logger.info("Splitting contigs")
        # TODO
        checkpoint[stage] = stage_params
        mark_checkpoint(checkpoint, checkpoints_path)
    else:
        logger.info("Contigs already split")

    # finding single-copy marker genes
    stage = "2_3"
    stage_params = [contigs]

    if True or should_run_step(stage_params, checkpoint, stage):
        logger.info("Searching for marker genes")
        # values in each entry of marker_contigs dict are must-not-links
        marker_gene_utils.scan_for_marker_genes(contigs, threads)
        marker_contigs, marker_contig_counts, contig_markers = \
            marker_gene_utils.get_contigs_with_marker_genes(contigs, 0.5, contig_length, 1000)
        checkpoint[stage] = stage_params
        mark_checkpoint(checkpoint, checkpoints_path)
    else:
        logger.info("Loading marker genes from previous computations")


    # compute k-mer vectors
    stage = "3_1"
    stage_params = [contigs, k_size]

    if should_run_step(stage_params, checkpoint, stage):
        logger.info("Computing k-mer vectors")
        run_kmers(reads_path, output, k_size, threads)
        checkpoint[stage] = stage_params
        mark_checkpoint(checkpoint, checkpoints_path)
        logger.info("Computing k-mer vectors complete")
    else:
        logger.info("K-mer vectors already computed")


    # compute coverage vectors
    stage = "3_2"
    stage_params = [contigs, bin_size, bin_count]

    if should_run_step(stage_params, checkpoint, stage):
        logger.info("Counting k-mers")
        # TODO
        checkpoint[stage] = stage_params
        mark_checkpoint(checkpoint, checkpoints_path)
        logger.info("Counting k-mers complete")
    else:
        logger.info("K-mer vectors already computed")


def run_reads_binning(args):
    # commong arguments for reads binning and contig binning
    reads_path = args.reads_path
    threads = args.threads
    bin_size = args.bin_size
    bin_count = args.bin_count
    k_size = args.k_size
    epochs = args.ae_epochs
    dims = args.ae_dims
    hidden = list(map(int, args.ae_hidden.split(",")))
    binreads = args.separate_reads
    cuda = args.cuda
    resume = args.resume
    min_cluster_size = max(args.min_bin_size, 1)
    iterations = max(args.bin_iterations, 1)
    output = args.output

    checkpoints_path = f"{output}/checkpoints"

    if not resume:
        checkpoint = {}
        mark_checkpoint(checkpoint, checkpoints_path)
    else:
        logger.info("Resuming the program from previous checkpoints")
        checkpoint = load_checkpoints(checkpoints_path)
        logger.debug(str(checkpoint))

    # computing k-mer vectors
    stage = "1_1"
    stage_params = [reads_path, k_size]

    if should_run_step(stage_params, checkpoint, stage):
        logger.info("Counting k-mers")
        run_kmers(reads_path, output, k_size, threads)
        checkpoint[stage] = stage_params
        mark_checkpoint(checkpoint, checkpoints_path)
        logger.info("Counting k-mers complete")
    else:
        logger.info("K-mer vectors already computed")

    # counting 15-mers
    stage = "1_2"
    stage_params = [reads_path]

    if should_run_step(stage_params, checkpoint, stage):
        logger.info("Counting 15-mers")
        run_15mer_counts(reads_path, output, threads)
        checkpoint[stage] = stage_params
        mark_checkpoint(checkpoint, checkpoints_path)
        logger.info("Counting 15-mers complete")
    else:
        logger.info("15-mers already counted")

    # computing coverage vectors
    stage = "2_1"
    stage_params = [reads_path, bin_size, bin_count]

    if should_run_step(stage_params, checkpoint, stage):
        logger.info("Computing 15-mer profiles")
        run_15mer_vecs(
            reads_path, output, bin_size, bin_count, threads)
        checkpoint[stage] = stage_params
        mark_checkpoint(checkpoint, checkpoints_path)
        logger.info("Computing 15-mer profiles complete")
    else:
        logger.info("Already computed 15-mer profiles complete")


    # numpy vectors
    stage = "3_1"
    stage_params = ['numpy']

    if should_run_step(stage_params, checkpoint, stage):
        logger.info("Profiles saving as numpy arrays")
        comp_profiles = np.array([np.array(list(map(float, line.strip().split()))) for line in open(
            f"{output}/profiles/com_profs") if len(line.strip()) > 0])
        cov_profiles = np.array([np.array(list(map(float, line.strip().split()))) for line in open(
            f"{output}/profiles/cov_profs") if len(line.strip()) > 0])

        np.save(f"{output}/profiles/com_profs", comp_profiles)
        np.save(f"{output}/profiles/cov_profs", cov_profiles)

        del comp_profiles
        del cov_profiles

        logger.info("Profiles saving as numpy arrays complete")
        checkpoint[stage] = stage_params
        mark_checkpoint(checkpoint, checkpoints_path)
    else:
        logger.info("Numpy arrays already computed")

    # VAE encode
    # TODO 
    constraints = None
    stage = "4_1"
    stage_params = [output,
                    dims,
                    hidden,
                    epochs,
                    constraints]

    if should_run_step(stage_params, checkpoint, stage):

        logger.info(f"VAE training information")
        logger.info(f"\tDimensions {dims}")
        logger.info(f"\tHidden Layers {hidden}")
        logger.info(f"\tEpochs {epochs}")

        ae_utils.vae_encode(
            output,
            dims,
            hidden,
            epochs,
            constraints,
            cuda)

        checkpoint[stage] = stage_params
        mark_checkpoint(checkpoint, checkpoints_path)
        logger.info(f"VAE training complete")
    else:
        logger.info(f"VAE already trained")

    # Must run content
    if constraints is None:
        cluster_utils.perform_binning(
            output, iterations, min_cluster_size, binreads, reads_path)
    else:
        cluster_utils.perform_binning_HDBSCAN(
            output, min_cluster_size, binreads, reads_path, threads)