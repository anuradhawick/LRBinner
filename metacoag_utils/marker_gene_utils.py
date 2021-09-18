#!/usr/bin/env python3

import os
import pathlib
import logging

# create logger
logger = logging.getLogger('LRBinner')


# Modified from SolidBin
def scan_for_marker_genes(contig_file, nthreads, hard=0):

    software_path = pathlib.Path(__file__).parent.absolute()

    fragScanURL = os.path.join(software_path.parent, 'auxiliary',
                               'FragGeneScan1.31', 'run_FragGeneScan.pl')
    hmmExeURL = os.path.join(software_path.parent, 'auxiliary', 'hmmer-3.3.2',
                             'src', 'hmmsearch')
    markerURL = os.path.join(software_path.parent, 'auxiliary', 'marker.hmm')

    logger.debug(markerURL)

    fragResultURL = contig_file+".frag.faa"
    hmmResultURL = contig_file+".hmmout"

    if not (os.path.exists(fragResultURL)):
        fragCmd = fragScanURL+" -genome="+contig_file+" -out="+contig_file + \
            ".frag -complete=0 -train=complete -thread="+str(nthreads)+" 1>" + \
            contig_file+".frag.out 2>"+contig_file+".frag.err"
        logger.debug("exec cmd: "+fragCmd)
        os.system(fragCmd)

    if os.path.exists(fragResultURL):
        if not (os.path.exists(hmmResultURL)):
            hmmCmd = hmmExeURL+" --domtblout "+hmmResultURL+" --cut_tc --cpu "+str(nthreads)+" " + \
                markerURL+" "+fragResultURL+" 1>"+hmmResultURL+".out 2>"+hmmResultURL+".err"
            logger.debug("exec cmd: "+hmmCmd)
            os.system(hmmCmd)

        else:
            logger.error("HMMER search failed! Path: " +
                         hmmResultURL + " does not exist.")
    else:
        logger.error("FragGeneScan failed! Path: " +
                     fragResultURL + " does not exist.")


# Get contigs containing marker genes
def get_contigs_with_marker_genes(
        contigs_file, mg_length_threshold, contig_lengths, min_length):

    marker_contigs = {}
    marker_contig_counts = {}
    contig_markers = {}

    with open(contigs_file+".hmmout", "r") as myfile:
        for line in myfile.readlines():
            if not line.startswith("#"):
                strings = line.strip().split()

                contig = strings[0]

                # Marker gene name
                marker_gene = strings[3]

                # Marker gene length
                marker_gene_length = int(strings[5])

                # Mapped marker gene length
                mapped_marker_length = int(strings[16]) - int(strings[15])

                # Get contig name
                name_strings = contig.split("_")
                name_strings = name_strings[:len(name_strings)-3]
                contig_name = "_".join(name_strings)

                contig_length = contig_lengths[contig_name]

                if contig_length >= min_length and mapped_marker_length > marker_gene_length*mg_length_threshold:

                    marker_repeated_in_contig = False

                    # Get marker genes in each contig
                    if contig_name not in contig_markers:
                        contig_markers[contig_name] = [marker_gene]
                    else:
                        if marker_gene not in contig_markers[contig_name]:
                            contig_markers[contig_name].append(marker_gene)

                    # Get contigs containing each marker gene
                    if marker_gene not in marker_contigs:
                        marker_contigs[marker_gene] = [contig_name]
                    else:
                        if contig_name not in marker_contigs[marker_gene]:
                            marker_contigs[marker_gene].append(contig_name)
                        else:
                            marker_repeated_in_contig = True

                    # Get contig counts for each marker
                    if marker_gene not in marker_contig_counts:
                        marker_contig_counts[marker_gene] = 1
                    else:
                        if not marker_repeated_in_contig:
                            marker_contig_counts[marker_gene] += 1

    return marker_contigs, marker_contig_counts, contig_markers


def count_contigs_with_marker_genes(marker_contig_counts):

    marker_frequencies = {}

    for marker in marker_contig_counts:

        if marker_contig_counts[marker] not in marker_frequencies:
            marker_frequencies[marker_contig_counts[marker]] = 1
        else:
            marker_frequencies[marker_contig_counts[marker]] += 1

    return marker_frequencies
