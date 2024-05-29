# Genetic and recombination maps

A genetic map for Bgt was obtained by Müller et al. 2019 [doi.org/10.1111/nph.15529](https://doi.org/10.1111/nph.15529)  by crossing two isolates and sequencing the genome of 117 progenies.

We retreived the file `GeneticMap_final_13052018.csv` from the [github page](https://github.com/MarionCMueller/AvrPm17/tree/main/Figure1) of the first author. This file containes a list of markers and their distance in cM, bu no information on their physical position. This information was provide by the authors of the paper above (file `Positions_THUN12X96224.csv`)

First we combined the genetic and physical information into one file:

```
# extract info from genetic maps 
head -n 1 ../data/GeneticMap_final_13052018.csv | sed 's/Isolate,//g' | sed 's/,/\n/g' > list_snp_ig_G
head -n 2 ../data/GeneticMap_final_13052018.csv| tail -n 1 | sed 's/^,1/1/g' | sed 's/,/\n/g'  > list_chr_G
head -n 3 ../data/GeneticMap_final_13052018.csv| tail -n 1 | sed 's/^,0/0/g' | sed 's/,/\n/g'  > list_rec_G

#extract physical positions to attach to genetic map
grep -w -f list_snp_ig_G ../data/Positions_THUN12X96224.csv  > Positions_THUN12X96224_present_in_G.csv

cut -d ";" -f 3 Positions_THUN12X96224_present_in_G.csv  > list_physical_position_G

paste  -d "," list_snp_ig_G list_chr_G list_rec_G list_physical_position_G  > THUN12x96224_genetic_map_in_cM_+_phy_distance
```

Second, we converted the genetic map in to a recombination map by calculating the local per base recombination rates:

```
python calc_rec_rate_per_bp.py -o THUN12x96224_bp_recombination_rates.txt -rec THUN12x96224_genetic_map_in_cM_+_phy_distance
```
The resulting file `THUN12x96224_bp_recombination_rates.txt` was then used in downstream analyses such as [fineStructure](../fineStructure/fineStructure.md) and [isoRelate](../isoRelate/isoRelate.md).

For imputation with Beagle we need a slightly different format

```
beagle code
```
And for Relate as well we need to modify the format and a separate file for each chromosome (available in this folder):

```
python calc_rec_rate_per_bp_Relate.py -rec THUN12x96224_genetic_map_in_cM_+_phy_distance -o THUN12x96224_4Relate
```
