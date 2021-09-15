import numpy as np
from collections import defaultdict
import torch
import matplotlib.pyplot as plt
import math
import random
import pickle
from tqdm import tqdm
import logging
import os
from Bio import SeqIO

logger = logging.getLogger('LRBinner')

# start code from vamb


def smaller_indices(distances, threshold):
    arr = distances.numpy()
    indices = (arr <= threshold).nonzero()[0]
    torch_indices = torch.from_numpy(indices)

    return torch_indices


def normalize(matrix, inplace=False):
    if isinstance(matrix, np.ndarray):
        matrix = torch.from_numpy(matrix)

    matrix = matrix.clone()

    # If any rows are kept all zeros, the distance function will return 0.5 to all points
    # inclusive itself, which can break the code in this module
    zeromask = matrix.sum(dim=1) == 0
    matrix[zeromask] = 1/matrix.shape[1]
    matrix /= (matrix.norm(dim=1).reshape(-1, 1) * (2 ** 0.5))
    return matrix


def calc_distances(matrix, index):
    "Return vecembeddedRoottor of cosine distances from rows of normalized matrix to given row."
    dists = 0.5 - matrix.matmul(matrix[index])
    dists[index] = 0.0  # avoid float rounding errors
    return dists


_DELTA_X = 0.005
_XMAX = 0.3

# This is the PDF of normal with µ=0, s=0.01 from -0.075 to 0.075 with intervals
# of DELTA_X, for a total of 31 values. We multiply by _DELTA_X so the density
# of one point sums to approximately one
_NORMALPDF = _DELTA_X * torch.Tensor(
    [2.43432053e-11, 9.13472041e-10, 2.66955661e-08, 6.07588285e-07,
     1.07697600e-05, 1.48671951e-04, 1.59837411e-03, 1.33830226e-02,
     8.72682695e-02, 4.43184841e-01, 1.75283005e+00, 5.39909665e+00,
     1.29517596e+01, 2.41970725e+01, 3.52065327e+01, 3.98942280e+01,
     3.52065327e+01, 2.41970725e+01, 1.29517596e+01, 5.39909665e+00,
     1.75283005e+00, 4.43184841e-01, 8.72682695e-02, 1.33830226e-02,
     1.59837411e-03, 1.48671951e-04, 1.07697600e-05, 6.07588285e-07,
     2.66955661e-08, 9.13472041e-10, 2.43432053e-11])


def calc_densities(histogram, cuda=False, pdf=_NORMALPDF):
    """Given an array of histogram, smoothes the histogram."""
    pdf_len = len(pdf)

    if cuda:
        histogram = histogram.cpu()

    densities = torch.zeros(len(histogram) + pdf_len - 1)
    for i in range(len(densities) - pdf_len + 1):
        densities[i:i+pdf_len] += pdf * histogram[i]

    densities = densities[15:-15]

    return densities

# end code from vamb


def find_valley_ratio(densities):
    peak_density = 0
    min_density = None
    peak_over = False
    success = False

    minima = None
    maxima = None
    early_minima = None

    x = 0
    for n, density in enumerate(densities):
        if not peak_over and density > peak_density:
            if x > 0.1:
                break
            peak_density = density
            maxima = x

        if not peak_over and density < peak_density:
            peak_over = True
            peak_density = density
            min_density = density
            minima = x

        if peak_over and density > min_density:
            break

        # Now find the minimum after the peak
        if peak_over and density < min_density:
            min_density = density
            minima = x
            if n != 0 and (densities[n-1]-densities[n])/(1/_DELTA_X) > 0.5:
                early_minima = x

            # break on platue
            if (densities[n-1]-densities[n])/(1/_DELTA_X) < 0.2:
                break

        x += _DELTA_X

    if not peak_over:
        return False, False, False, False

    if early_minima is None:
        early_minima = minima

    return min_density/peak_density, maxima, early_minima, minima


def get_cluster_center(matrix, seed):
    distances = calc_distances(matrix, seed)
    histogram = torch.histc(distances, math.ceil(_XMAX/_DELTA_X), 0, _XMAX)
    histogram[0] -= 1
    densities = calc_densities(histogram)
    ratio, chosen_peak, chosen_minima, chosen_tail = find_valley_ratio(
        densities)

    # To make my life easy debugging I'm gonna keep this!
    # import matplotlib.pyplot as plt
    # print(distances.sort())
    # plt.figure()
    # plt.plot(np.arange(0, 2, _DELTA_X), calc_densities(torch.histc(distances, math.ceil(2/_DELTA_X), 0, 2)))
    # plt.show()
    # plt.close()

    if not chosen_peak or ratio > 0.5:
        return False, False, False, False, False

    from_x, to_x = chosen_peak-_DELTA_X*5, chosen_peak+_DELTA_X*5
    chosen_points = [p for p, x in enumerate(
        distances.numpy()) if to_x > x > from_x]

    if len(chosen_points) < 100:
        return False, False, False, False, False

    sample_size = min(1000, max(100, len(chosen_points) * 0.01))
    sample_size = int(sample_size)
    sampled_points = random.sample(chosen_points, sample_size)

    ratio = 10000
    best_point = None
    best_densities = None
    distance_cache = None
    tail = None
    minima = None
    maxima = None

    for p in sampled_points:
        distances = calc_distances(matrix, p)
        histogram = torch.histc(distances, math.ceil(_XMAX/_DELTA_X), 0, _XMAX)
        histogram[0] -= 1
        densities = calc_densities(histogram)

        new_ratio, new_maxima, new_minima, new_tail = find_valley_ratio(
            densities)

        if new_ratio and new_ratio < ratio:
            ratio = new_ratio
            best_point = p
            best_densities = densities
            distance_cache = distances
            tail = new_tail
            minima = new_minima
            maxima = new_maxima

    return best_point, distance_cache, maxima, minima, tail


def cluster_points(latent, iterations, min_cluster_size):
    matrix = normalize(latent)
    clusters = defaultdict(list)
    read_ids = np.arange(len(matrix))
    read_ids_ref = np.arange(len(matrix))

    if iterations != 0:
        for x in tqdm(range(iterations), total=iterations, desc="Performing iterations"):
            if len(read_ids) < min_cluster_size * 0.6:
                break
            random_point = random.choice(read_ids)
            best_point, distance_cache, maxima, minima, tail = get_cluster_center(
                matrix, random_point)

            if tail:
                cluster_pts = smaller_indices(distance_cache, tail)
                removables = smaller_indices(distance_cache, tail)
                removables_idx = set(read_ids_ref[removables])
                clusters[x] = set(read_ids_ref[cluster_pts])

                new_read_ids_ref = np.array(
                    [y for y in read_ids_ref if y not in removables_idx])
                new_matrix = np.delete(matrix.numpy(), removables, axis=0)
                new_read_ids = np.arange(len(new_read_ids_ref))
                matrix = torch.from_numpy(new_matrix)
                read_ids = new_read_ids
                read_ids_ref = new_read_ids_ref
    else:
        x = 0

        while True:
            if len(read_ids) < min_cluster_size * 0.1:
                break

            finish_search = True
            random_candidates = list(read_ids)
            random.shuffle(random_candidates)

            for random_point in tqdm(random_candidates, desc="Performing exhaustive search!"):
                best_point, distance_cache, maxima, minima, tail = get_cluster_center(
                    matrix, random_point)

                if tail:
                    cluster_pts = smaller_indices(distance_cache, tail)
                    removables = smaller_indices(distance_cache, tail)
                    removables_idx = set(read_ids_ref[removables])
                    clusters[x] = set(read_ids_ref[cluster_pts])
                    x += 1

                    new_read_ids_ref = np.array(
                        [y for y in read_ids_ref if y not in removables_idx])
                    new_matrix = np.delete(matrix.numpy(), removables, axis=0)
                    new_read_ids = np.arange(len(new_read_ids_ref))
                    matrix = torch.from_numpy(new_matrix)
                    read_ids = new_read_ids
                    read_ids_ref = new_read_ids_ref
                    finish_search = False
                    break
            if finish_search:
                break

    return clusters


def normal(val, mean, std):
    a = np.sqrt(2*np.pi) * std
    b = -0.5 * np.square((val-mean)/std)
    b = np.exp(b)
    c = b/a + 0.0000001
    pdf = np.sum(np.log(c))

    return pdf


def perform_binning(output, iterations, min_cluster_size, binreads, reads):
    latent = np.load(f'{output}/latent.npy')
    logger.info("Clustering algorithm running")
    clusters = cluster_points(latent, iterations, min_cluster_size)
    clusters_output = {}
    logger.info(f"Detected {len(clusters)}")

    for k, v in clusters.items():
        if len(v) > min_cluster_size:
            clusters_output[len(clusters_output)] = list(map(int, v))
    logger.info(
        f"Detected {len(clusters_output)} clusters with more than {min_cluster_size} points")
    cluster_profiles = {}
    classified_reads = []

    logger.info("Building profiles")

    comp_profiles = np.load(f"{output}/profiles/com_profs.npy")
    cov_profiles = np.load(f"{output}/profiles/cov_profs.npy")

    for k, rs in clusters_output.items():
        vecs = []
        for r in rs:
            classified_reads.append(r)
            vecs.append(np.concatenate(
                [comp_profiles[r], cov_profiles[r]], axis=0))
        vecs = np.array(vecs)
        cluster_profiles[k] = {
            'mean': vecs.mean(axis=0),
            'std': vecs.std(axis=0),
        }
    classified_reads = set(classified_reads)

    all_reads = set(range(len(comp_profiles)))
    unclassified_reads = all_reads - classified_reads

    logger.debug(f"Unclassified points to cluster {len(unclassified_reads)}")
    logger.info(f"Binning unclassified reads")
    for r in tqdm(unclassified_reads):
        max_p = float('-inf')
        best_c = None

        for k, v in cluster_profiles.items():
            p = normal(np.concatenate(
                [comp_profiles[r], cov_profiles[r]], axis=0), v['mean'], v['std'])

            if p > max_p:
                max_p = p
                best_c = k

        if best_c is not None:
            clusters_output[best_c].append(r)

    logger.info(f"Binning complete with {len(clusters_output)} bins")
    pickle.dump(clusters_output, open(f"{output}/binning_result.pkl", "wb+"))

    # separating reads into bins
    if binreads:
        bin_files = {}
    read_bin = {}

    if binreads:
        if os.path.isdir(f"{output}/binned_reads"):
            os.rmdir(f"{output}/binned_reads")

        if not os.path.isdir(f"{output}/binned_reads"):
            os.makedirs(f"{output}/binned_reads")

    for k, v in clusters_output.items():
        if binreads:
            if not os.path.isdir(f"{output}/bin-{k}/"):
                os.makedirs(f"{output}/bin-{k}/")

            bin_files[k] = open(f"{output}/bin-{k}/reads.fasta", "w+")

        for r in v:
            read_bin[r] = k

    binout = open(f"{output}/bins.txt", "w+")
    lenout = open(f"{output}/lengths.txt", "w+")
    fmt = "fasta" if reads.split(
        '.')[-1].lower() in ["fasta", "fna", "fa"] else "fastq"

    for r, record in enumerate(SeqIO.parse(reads, fmt)):
        binout.write(f"{read_bin[r]}\n")
        lenout.write(f"{len(record.seq)}\n")
        if binreads:
            bin_files[read_bin[r]].write(f">read-{r}\n")
            bin_files[read_bin[r]].write(f"{record.seq}\n")

    binout.close()
    lenout.close()

    if binreads:
        for k, f in bin_files.items():
            f.close()

def perform_binning_HDBSCAN(output, min_cluster_size, binreads, reads):
    import umap
    from hdbscan import HDBSCAN
    
    latent = np.load(f'{output}/latent.npy')
    logger.info("Clustering using HDBSCAN running")

    sidx = random.sample(range(len(latent)), 50000)
    sampled_data = latent[sidx]

    embds = umap.UMAP().fit_transform(sampled_data)
    labels = HDBSCAN(min_cluster_size=500).fit_predict(embds)

    clusters = defaultdict(list)

    for i, c in zip(sidx, labels):
        if c != -1:
            clusters[c].append(i)

    logger.info(f"Detected {len(clusters)}")

    clusters_output = {}

    for k, v in clusters.items():
        if len(v) > 50000 * min_cluster_size/len(latent):
            clusters_output[len(clusters_output)] = list(map(int, v))
    logger.info(
        f"Detected {len(clusters_output)} clusters with more than {min_cluster_size} points")
    cluster_profiles = {}
    classified_reads = []

    logger.info("Building profiles")

    comp_profiles = np.load(f"{output}/profiles/com_profs.npy")
    cov_profiles = np.load(f"{output}/profiles/cov_profs.npy")

    for k, rs in clusters_output.items():
        vecs = []
        for r in rs:
            classified_reads.append(r)
            vecs.append(np.concatenate(
                [comp_profiles[r], cov_profiles[r]], axis=0))
        vecs = np.array(vecs)
        cluster_profiles[k] = {
            'mean': vecs.mean(axis=0),
            'std': vecs.std(axis=0),
        }
    classified_reads = set(classified_reads)

    all_reads = set(range(len(comp_profiles)))
    unclassified_reads = all_reads - classified_reads

    logger.debug(f"Unclassified points to cluster {len(unclassified_reads)}")
    logger.info(f"Binning unclassified reads")
    for r in tqdm(unclassified_reads):
        max_p = float('-inf')
        best_c = None

        for k, v in cluster_profiles.items():
            p = normal(np.concatenate(
                [comp_profiles[r], cov_profiles[r]], axis=0), v['mean'], v['std'])

            if p > max_p:
                max_p = p
                best_c = k

        if best_c is not None:
            clusters_output[best_c].append(r)

    logger.info(f"Binning complete with {len(clusters_output)} bins")
    pickle.dump(clusters_output, open(f"{output}/binning_result.pkl", "wb+"))

    # separating reads into bins
    if binreads:
        bin_files = {}
    read_bin = {}

    if binreads:
        if os.path.isdir(f"{output}/binned_reads"):
            os.rmdir(f"{output}/binned_reads")

        if not os.path.isdir(f"{output}/binned_reads"):
            os.makedirs(f"{output}/binned_reads")

    for k, v in clusters_output.items():
        if binreads:
            if not os.path.isdir(f"{output}/bin-{k}/"):
                os.makedirs(f"{output}/bin-{k}/")

            bin_files[k] = open(f"{output}/bin-{k}/reads.fasta", "w+")

        for r in v:
            read_bin[r] = k

    binout = open(f"{output}/bins.txt", "w+")
    lenout = open(f"{output}/lengths.txt", "w+")
    fmt = "fasta" if reads.split(
        '.')[-1].lower() in ["fasta", "fna", "fa"] else "fastq"

    for r, record in enumerate(SeqIO.parse(reads, fmt)):
        binout.write(f"{read_bin[r]}\n")
        lenout.write(f"{len(record.seq)}\n")
        if binreads:
            bin_files[read_bin[r]].write(f">read-{r}\n")
            bin_files[read_bin[r]].write(f"{record.seq}\n")

    binout.close()
    lenout.close()

    if binreads:
        for k, f in bin_files.items():
            f.close()