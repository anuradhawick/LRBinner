import pickle
import os
import logging
import sys
from collections import defaultdict
from Bio import SeqIO
import numpy as np
import random

from mbcclr_utils.runners_utils import *
from mbcclr_utils import ae_utils
from mbcclr_utils import cluster_utils

logger = logging.getLogger('LRBinner')

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
    separate = args.separate
    cuda = args.cuda
    resume = args.resume
    contigs = args.contigs
    output = args.output
    cuda = args.cuda

    checkpoints_path = f"{output}/checkpoints"

    if not resume:
        checkpoint = Checkpointer(checkpoints_path)
    else:
        logger.info("Resuming the program from previous checkpoints")
        checkpoint = Checkpointer(checkpoints_path, True)
        logger.debug(checkpoint)

    stage = "1_1"
    stage_params = ['contigs_binning']

    if checkpoint.should_run_step(stage, stage_params):
        checkpoint.log(stage, stage_params)

    # compute contig lengths
    stage = "2_1"
    stage_params = [contigs]

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Computing contig lengths")
        contig_length = {}
        contig_id_idx = {}
        contig_idx_id = {}

        for record in SeqIO.parse(contigs, "fasta"):
            contig_length[record.id] = len(record.seq)
            contig_idx_id[len(contig_id_idx)] = record.id
            contig_id_idx[record.id] = len(contig_id_idx)
        
        pickle.dump(contig_length, open(f"{output}/profiles/contig_lengths.pkl", "wb+"))
        pickle.dump(contig_id_idx, open(f"{output}/profiles/contig_id_idx.pkl", "wb+"))
        pickle.dump(contig_idx_id, open(f"{output}/profiles/contig_idx_id.pkl", "wb+"))
        checkpoint.log(stage, stage_params)
    else:
        logger.info("Loading contig lengths")
        contig_length = pickle.load(open(f"{output}/profiles/contig_lengths.pkl", "rb"))
        contig_id_idx = pickle.load(open(f"{output}/profiles/contig_id_idx.pkl", "rb"))
        contig_idx_id = pickle.load(open(f"{output}/profiles/contig_idx_id.pkl", "rb"))

    # finding single-copy marker genes
    stage = "2_2"
    stage_params = [contigs]

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Searching for marker genes")
        # values in each entry of marker_contigs dict are must-not-links
        marker_gene_utils.scan_for_marker_genes(contigs, output, threads)
        marker_contigs, _, _ = \
            marker_gene_utils.get_contigs_with_marker_genes(output, 0.5, contig_length, 1000)
        pickle.dump(marker_contigs, open(f"{output}/profiles/marker_contigs.pkl", "wb+"))
        
        checkpoint.log(stage, stage_params)
        logger.info("Searching for marker genes complete")
    else:
        marker_contigs = pickle.load(open(f"{output}/profiles/marker_contigs.pkl", "rb"))
        logger.info("Loading marker genes from previous computations")

    # splitting contigs
    stage = "2_3"
    stage_params = [contigs]

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Splitting contigs")
        contig_groups, fragment_parent = split_contigs(contigs, output)
        must_link_pairs = []
        must_not_link_contig_pairs = []
        must_not_link_pairs = []

        # for _, frag_id in contig_groups.items():
        #     for i in random.sample(frag_id, min(10, len(frag_id))):
        #         for j in random.sample(frag_id, min(10, len(frag_id))):
        #             if j!=i:
        #                 must_link_pairs.append([i, j])
        
        for _, contig_ids in marker_contigs.items():
            contig_ids = sorted([contig_id_idx[x] for x in contig_ids]) 

            for c1 in contig_ids:
                for c2 in contig_ids:
                    if c2>=c1:
                        break
                    must_not_link_contig_pairs.append([c1, c2])
        
        for c1idx, c2idx in must_not_link_contig_pairs:
            g1 = contig_groups[contig_idx_id[c1idx]]
            g2 = contig_groups[contig_idx_id[c2idx]]

            g11 = random.sample(g1, min(10, len(g1)))
            g22 = random.sample(g2, min(10, len(g2)))
            
            for gg1 in g11:
                for gg2 in g22:
                    must_not_link_pairs.append([gg1, gg2])

        pickle.dump(must_link_pairs, open(f"{output}/profiles/must_link_pairs.pkl", "wb+"))
        pickle.dump(must_not_link_pairs, open(f"{output}/profiles/must_not_link_pairs.pkl", "wb+"))
        pickle.dump(contig_groups, open(f"{output}/profiles/contig_groups.pkl", "wb+"))
        pickle.dump(fragment_parent, open(f"{output}/profiles/fragment_parent.pkl", "wb+"))
        
        checkpoint.log(stage, stage_params)
        logger.info("Splitting contigs completed")
    else:
        must_link_pairs = pickle.load(open(f"{output}/profiles/must_link_pairs.pkl", "rb"))
        must_not_link_pairs = pickle.load(open(f"{output}/profiles/must_not_link_pairs.pkl", "rb"))
        contig_groups = pickle.load(open(f"{output}/profiles/contig_groups.pkl", "rb"))
        fragment_parent = pickle.load(open(f"{output}/profiles/fragment_parent.pkl", "rb"))
        logger.info("Contigs already split")

    # counting 15-mers
    stage = "2_4"
    stage_params = [reads_path]

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Counting 15-mers")
        run_15mer_counts(reads_path, output, threads)
        
        checkpoint.log(stage, stage_params)
        logger.info("Counting 15-mers complete")
    else:
        logger.info("15-mer counting already performed")

    # compute k-mer vectors
    stage = "3_1"
    fragments_path = f"{output}/fragments/contigs.fasta"
    stage_params = [fragments_path, k_size]

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Computing k-mer vectors")
        run_kmers(fragments_path, output, k_size, threads)
        checkpoint.log(stage, stage_params)
        logger.info("Computing k-mer vectors complete")
    else:
        logger.info("K-mer vectors already computed")

    # compute coverage vectors
    stage = "4_1"
    stage_params = [fragments_path, bin_size, bin_count]

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Generating coverage vectors")
        run_15mer_vecs(fragments_path, output, bin_size, bin_count, threads)
        checkpoint.log(stage, stage_params)
        logger.info("Generating coverage vectors complete")
    else:
        logger.info("Coverage vectors already computed")

    # numpy vectors
    stage = "5_1"
    stage_params = ['numpy']

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Profiles saving as numpy arrays")
        comp_profiles = np.array([np.array(list(map(float, line.strip().split()))) for line in open(
            f"{output}/profiles/com_profs") if len(line.strip()) > 0])
        cov_profiles = np.array([np.array(list(map(float, line.strip().split()))) for line in open(
            f"{output}/profiles/cov_profs") if len(line.strip()) > 0])

        np.save(f"{output}/profiles/com_profs", comp_profiles)
        np.save(f"{output}/profiles/cov_profs", cov_profiles)
        data_length = len(cov_profiles)
        del comp_profiles
        del cov_profiles

        logger.info("Profiles saving as numpy arrays complete")
        
        checkpoint.log(stage, stage_params)
    else:
        data_length = len(np.load(f"{output}/profiles/com_profs.npy"))
        logger.info("Numpy arrays already computed")

    # VAE encode
    stage = "6_1"
    stage_params = [output,
                    dims,
                    hidden,
                    epochs,
                    len(must_link_pairs),
                    len(must_not_link_pairs)]

    if checkpoint.should_run_step(stage, stage_params):

        logger.info(f"VAE training information")
        logger.info(f"\tDimensions {dims}")
        logger.info(f"\tHidden Layers {hidden}")
        logger.info(f"\tEpochs {epochs}")

        # constraints = None
        constraints = {'ml': must_link_pairs, 
                       'mnl': must_not_link_pairs,
                       'size': data_length}
        
        logger.info(f"Contig split must link pairs   {len(must_link_pairs):10}")
        logger.info(f"Single copy marker genes pairs {len(must_not_link_pairs):10}")

        ae_utils.vae_encode(
            output,
            dims,
            hidden,
            epochs,
            constraints,
            cuda)

        checkpoint.log(stage, stage_params)
        logger.info(f"VAE training complete")
    else:
        logger.info(f"VAE already trained")

    # Must run content
    cluster_utils.perform_contig_binning_HDBSCAN(
            output, fragment_parent, separate, contigs, threads)

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
    separate = args.separate
    cuda = args.cuda
    resume = args.resume
    min_cluster_size = max(args.min_bin_size, 1)
    iterations = max(args.bin_iterations, 0)
    output = args.output

    checkpoints_path = f"{output}/checkpoints"

    if not resume:
        checkpoint = Checkpointer(checkpoints_path)
    else:
        logger.info("Resuming the program from previous checkpoints")
        checkpoint = Checkpointer(checkpoints_path, True)
        logger.debug(checkpoint)

    # computing k-mer vectors
    stage = "1_1"
    stage_params = [reads_path, k_size]

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Counting k-mers")
        run_kmers(reads_path, output, k_size, threads)
        
        checkpoint.log(stage, stage_params)
        logger.info("Counting k-mers complete")
    else:
        logger.info("K-mer vectors already computed")

    # counting 15-mers
    stage = "1_2"
    stage_params = [reads_path]

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Counting 15-mers")
        run_15mer_counts(reads_path, output, threads)
        
        checkpoint.log(stage, stage_params)
        logger.info("Counting 15-mers complete")
    else:
        logger.info("15-mers already counted")

    # computing coverage vectors
    stage = "2_1"
    stage_params = [reads_path, bin_size, bin_count]

    if checkpoint.should_run_step(stage, stage_params):
        logger.info("Computing 15-mer profiles")
        run_15mer_vecs(
            reads_path, output, bin_size, bin_count, threads)
        
        checkpoint.log(stage, stage_params)
        logger.info("Computing 15-mer profiles complete")
    else:
        logger.info("Already computed 15-mer profiles complete")


    # numpy vectors
    stage = "3_1"
    stage_params = ['numpy']

    if checkpoint.should_run_step(stage, stage_params):
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
        
        checkpoint.log(stage, stage_params)
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

    if checkpoint.should_run_step(stage, stage_params):

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

        checkpoint.log(stage, stage_params)
        logger.info(f"VAE training complete")
    else:
        logger.info(f"VAE already trained")

    # Must run content
    if constraints is None:
        cluster_utils.perform_binning(
            output, iterations, min_cluster_size, separate, reads_path)
    else:
        cluster_utils.perform_binning_HDBSCAN(
            output, min_cluster_size, separate, reads_path, threads)