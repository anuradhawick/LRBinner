from tabulate import tabulate
import argparse

parser = argparse.ArgumentParser(description='Separate reads in to bins.')

parser.add_argument(
    '--truth', '-t', help="Path to text file with grounds truth. Items with label Unknown will be ignored.", type=str, required=True)
parser.add_argument(
    '--bins', '-b', help="Path of bins.txt file from LRBinner.", type=str, required=True)

args = parser.parse_args()

truth = args.truth
bins = args.bins


def clusters_table(clusters, truth):
    clx = set(clusters)
    trx = set(truth)
    c_map = {k: v for v, k in enumerate(clx)}
    t_map = {k: v for v, k in enumerate(trx)}

    matrix = [["_"] + [f"Bin-{c_map[x]}_({x})" for x in list(clx)]]
    matrix += [[x] + [0 for i in range(len(clx))] for x in trx]

    mat = [[0 for i in range(len(clx))] for x in trx]

    for c, t in zip(clusters, truth):
        matrix[t_map[t] + 1][c_map[c] + 1] += 1
        mat[t_map[t]][c_map[c]] += 1

    matT = [[mat[j][i] for j in range(len(mat))] for i in range(len(mat[0]))]

    recall = sum([max(row) for row in mat])/sum([sum(row) for row in mat])
    precision = sum([max(row) for row in matT])/sum([sum(row) for row in matT])

    print(tabulate(matrix, tablefmt="plain"))
    print()
    print(f"Precision\t{recall*100:10.2f}")
    print(f"Recall    \t{precision*100:10.2f}")
    print(f"F1-Score  \t{(2 * recall*precision/(recall+precision))*100:10.2f}")


truth = open(truth).read().strip().split("\n")
bins = open(bins).read().strip().split("\n")

bins_f = [b for b, t in zip(bins, truth) if t.lower() != "unknown"]
truth_f = [t for b, t in zip(bins, truth) if t.lower() != "unknown"]

clusters_table(bins_f, truth_f)
