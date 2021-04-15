# cgMLST bootstrap trees
This is a script to create bootstrap trees for categorical data (like core genome MLST) from [Ridom SeqSphere](https://www.ridom.de/seqsphere/)

NOTE: currently the categorical data is required to be integer interpretable, all other values and columns are ignored for distance caluclation.
Sample ID column is identified if exported from SeqSphere (name: "Sample ID" - not case sensitive), otherwise column 0 is taken as ID

# requirements
- biopython
- pandas
- numpy

# usage

```
usage: bootstrap_alleles.py [-h] [-i TABLE] [-o OUTPUT] [-t TIMES]
                            [--tree TREE] [--upgma] [-d DELIM]

Routine to create bootstrap trees from categorical data likecore genome MLST

optional arguments:
  -h, --help            show this help message and exit
  -i TABLE, --input_table TABLE
                        Input table as .csv e.g. exported from SeqSphere
  -o OUTPUT, --output OUTPUT
                        destination to write newick tree [./tree.nwk]
  -t TIMES, --times TIMES
                        Number of bootstrap samples taken [100]
  --tree TREE           Newick tree file to add support values, if not
                        provided a new one is created
  --upgma               create UPGMA instead of NJ trees
  -d DELIM, --delim DELIM
                        Field delimiter [;]
```
