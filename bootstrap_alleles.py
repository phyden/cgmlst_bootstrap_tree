#!/usr/bin/env python
# coding: utf-8

import argparse
import numpy as np
import pandas as pd
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import Phylo
from Bio.Phylo.Consensus import get_support

parser = argparse.ArgumentParser(description="Routine to create bootstrap trees from categorical data like"
                                             " core genome MLST")
parser.add_argument("-i", "--input_table", dest="table",
                    help="Input table as .csv e.g. exported from SeqSphere")
parser.add_argument("-o", "--output", dest="output", default="tree.nwk",
                    help="destination to write newick tree [./tree.nwk]")
parser.add_argument("-t", "--times", dest="times", type=int, default=100,
                    help="Number of bootstrap samples taken [100]")
parser.add_argument("--tree", dest="tree", default=None,
                    help="Newick tree file to add support values, if not provided a new one is created")
parser.add_argument("--upgma", action="store_true", default=False,
                    help="create UPGMA instead of NJ trees")
parser.add_argument("-d", "--delim", dest="delim", type=str, default=";",
                    help="Field delimiter [;]")


def get_distance_matrix(array, names):
    d2 = array[:, np.newaxis, :]
    d2_t = d2.reshape((d2.shape[1], d2.shape[0], d2.shape[2]))
    dist = 1 - np.isclose(d2, d2_t)
    return to_bio_matrix(np.mean(dist, axis=2), names)


def to_bio_matrix(array, names):
    matrix = []
    for i in range(0, array.shape[0]):
        matrix.append(array[i, :].tolist()[:i + 1])
    return DistanceMatrix(names, matrix)


def bootstrap(array, times):
    i = 0
    while i < times:
        i += 1
        yield array[:, np.random.choice(array.shape[1], int(array.shape[1] * 0.9), replace=False)]


def main():
    args = parser.parse_args()

    table = args.table
    if not table:
        print("No input table provided. This argument is required\n")
        parser.print_help()
        exit(0)

    df = pd.read_csv(table, sep=args.delim, quotechar='"', index_col=False)
    for c in df.columns:
        if c.lower() == "sample id":
            df.set_index(c)
        else:
            df[c] = pd.to_numeric(df[c], errors="coerce", downcast="integer")

    df.dropna(axis=1, how="all", inplace=True)
    data = df.values.astype(int)
    names = df.index.tolist()

    constructor = DistanceTreeConstructor()
    if args.upgma:
        const = constructor.upgma
    else:
        const = constructor.nj

    if args.tree:
        main_tree = Phylo.read(args.tree, "newick")
    else:
        main_tree = const(get_distance_matrix(data, names))

    trees = [const(get_distance_matrix(sub, names)) for sub in bootstrap(data, args.times)]
    supp_tree = get_support(main_tree, trees)

    # remove labels of inner clades
    for clade in supp_tree.get_nonterminals():
        clade.name = ""

    Phylo.write(supp_tree, args.output, "newick")


if __name__ == '__main__':
    main()
