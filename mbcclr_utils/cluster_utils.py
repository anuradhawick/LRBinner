import numpy as np
from collections import defaultdict
import random
import seaborn as sns
import torch
import matplotlib.pyplot as plt 
import math
import random
from scipy.signal import find_peaks
import pickle
from tqdm import tqdm
from sklearn.metrics.cluster import adjusted_rand_score


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
    dists[index] = 0.0 # avoid float rounding errors
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
            if n!=0 and (densities[n-1]-densities[n])/(1/_DELTA_X) > 0.5:
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

def get_cluster_center(matrix, seed, peak_sampling_size):
    distances = calc_distances(matrix, seed)
    histogram = torch.histc(distances, math.ceil(_XMAX/_DELTA_X), 0, _XMAX)
    histogram[0] -= 1 
    densities = calc_densities(histogram)
    ratio, chosen_peak, chosen_minima, chosen_tail = find_valley_ratio(densities)
        
    if not chosen_peak or ratio > 0.5:
        return False, False, False, False, False

    from_x, to_x = chosen_peak-_DELTA_X*5, chosen_peak+_DELTA_X*5
    chosen_points = [p for p, x in enumerate(distances.numpy()) if to_x>x>from_x]

    if len(chosen_points) < peak_sampling_size:
        return False, False, False, False, False
    
    sampled_points = random.sample(chosen_points, min_cluster_size)
    
    
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
        
        new_ratio, new_maxima, new_minima, new_tail = find_valley_ratio(densities)  
        
        if new_ratio and new_ratio < ratio:
            ratio = new_ratio
            best_point = p
            best_densities = densities  
            distance_cache = distances
            tail = new_tail
            minima = new_minima
            maxima = new_maxima

    return best_point, distance_cache, maxima, minima, tail    

def cluster_points(latent, iterations, min_cluster_size=10000, peak_sampling_size=1000):
    matrix = normalize(latent)
    clusters = defaultdict(list)
    read_ids = np.arange(len(matrix))
    read_ids_ref = np.arange(len(matrix))

    for x in range(iterations):
        if len(read_ids) < min_cluster_size:
            break
        random_point = random.choice(read_ids)
        best_point, distance_cache, maxima, minima, tail = get_cluster_center(matrix, random_point, peak_sampling_size)
        

        if tail:
            cluster_pts = smaller_indices(distance_cache, tail)
            
            removables = smaller_indices(distance_cache, tail)
            removables_idx = set(read_ids_ref[removables])       
            clusters[x] = set(read_ids_ref[cluster_pts])
            
            new_read_ids_ref = np.array([y for y in read_ids_ref if y not in removables_idx])
            new_matrix = np.delete(matrix.numpy(), removables, axis=0)
            
            new_read_ids = np.arange(len(new_read_ids_ref))

            matrix = torch.from_numpy(new_matrix)
            read_ids = new_read_ids
            read_ids_ref = new_read_ids_ref
      
    return clusters

def normal(val, mean, std):
    a = np.sqrt(2*np.pi) * std
    b = -0.5 * np.square((val-mean)/std)
    b = np.exp(b)
    c = b/a + 0.0000001
    pdf = np.sum(np.log(c))

    return pdf

def perform_binning(save_path, latent, iterations, min_cluster_size=10000, peak_sampling_size=1000):
    clusters = cluster_points(latent, iterations, min_cluster_size/2, peak_sampling_size=1000)
    clusters_output = {}

    for k, v in clusters.items():
        if len(v) > min_cluster_size:
            clusters_output[len(clusters_output)] = list(map(int, v))

    cluster_profiles = {}
    classified_reads = []

    for k, rs in clusters_output.items():
        vecs = []
        for r in rs:
            classified_reads.append(r)
            vecs.append(np.concatenate([comp_profiles[r], cov_profiles[r]], axis=0))
        vecs = np.array(vecs)
        cluster_profiles[k] = {
            'mean': vecs.mean(axis=0),
            'std': vecs.std(axis=0),
        }
    classified_reads = set(classified_reads)

    all_reads = set(range(len(ground_truth)))
    unclassified_reads = all_reads - classified_reads

    for r in tqdm(unclassified_reads):
        max_p = float('-inf')
        best_c = None
        
        for k, v in cluster_profiles.items():
            p = normal(np.concatenate([comp_profiles[r], cov_profiles[r]], axis=0), v['mean'], v['std'])
            
            if p > max_p:
                max_p = p
                best_c = k
        
        if best_c is not None:
            clusters_output[best_c].append(r)

    pickle.dump(clusters_output, open(f"{save_path}/binning_result.pkl", "wb+"))