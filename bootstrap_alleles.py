#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import Phylo
from Bio.Phylo.Consensus import get_support
import matplotlib.pyplot as plt

# input table:
# should be col 1: sample id, rest integer (allele numbers)
table="./Human_oesterreich-2021.csv"

def get_distance_matrix(array, names):
    d2 = array[:, np.newaxis,:]
    d2_t = d2.reshape((d2.shape[1],d2.shape[0],d2.shape[2]))
    dist = 1 - np.isclose(d2, d2_t)
    return to_bio_matrix(np.mean(dist, axis=2), names)


def to_bio_matrix(array, names):
    matrix = []
    for i in range(0, array.shape[0]):
        matrix.append(array[i,:].tolist()[:i+1])
    return DistanceMatrix(names, matrix)
    

def bootstrap(array, times):
    i = 0
    while i < times:
        i += 1
        yield array[:, np.random.choice(array.shape[1], int(array.shape[1] // ( 1 / 0.9 )), replace=False)]


df = pd.read_csv(table, sep=';', quotechar='"', index_col=0)
for c in df.columns:
    df[c] = pd.to_numeric(df[c], errors="coerce", downcast="integer")

data = df.values.astype(int)
names = df.index.tolist()

constructor = DistanceTreeConstructor()
main_tree = constructor.nj(m)

main_tree = constructor.nj(get_distance_matrix(data, names))


trees = [constructor.nj(get_distance_matrix(sub, names)) for sub in bootstrap(data, 100)]

supp_tree = get_support(main_tree, trees)

Phylo.draw(supp_tree)
plt.show()

Phylo.write(supp_tree, "mytree.nwk", "newick")