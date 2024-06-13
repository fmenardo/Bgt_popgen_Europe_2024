# Sparg

We used sparg ([Osmond and Coop 2021](https://doi.org/10.1038/s41588-019-0484-x)) to estimate dispersal and infer ancestor locations for the dataset [Europe+_recent](../Datasets/Datasets.md).At the moment of writing the manuscript describing sparg was undergoing revisions. Therefore we include here the relevant sparg code obtained from the [sparg repository](https://github.com/mmosmond/sparg-ms/tree/revision) 

We divided the sparg analysis in 4 steps. 
With the first step we parsed the tree sequences output of Relate to select every 250th tree in the tree sequence and resample its branch lengths (and therefore coalescent times) 1000 times. All the input files were produced earlier with [Relate](../Relate/Relate.md).

```
mkdir data_random
python sparg_step1.py -in ../Relate/2022_2023_theta_popsize_chr$CHROMOSOME -o data_random/random250 -coal ../Relate/2022_2023_theta_popsize.coal -chr $CHROMOSOME -mut 0.0000005 -ns 1000 -min_snps 100 -min_win 50000 -nth_tree

```
In the second step we use sparg to process shared times and coalescent times, using a time cutoff of 100 generations.

```
python sparg_step2.py -coal ../Relate/2022_2023_theta_popsize.coal -tcut 100 -data_d data_random/ -c 10

```
The third step calculate the composite dispersal rate

```
python sparg_step3.py -meta 2022+before2022+2023+ncsu_metadata+fs+admxK9_03052024.csv -samples ../Relate/2022_2023.sample -ns 1000 -out random_250 -data_d data_random/
```

The fourth step we estimated abcestor locations for all samples at 1,2,3,4,5,10,25,50 and 100 generation before present.

```
python sparg_step4.py -tcut 100 -mle random_250.1000_trees.tcut.mle.npy -ns 1000 -data_d data_random/ -c 16
```
