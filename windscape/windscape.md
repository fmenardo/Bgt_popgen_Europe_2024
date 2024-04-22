# Windscape
First the wind data was downloaded from the [Climate System Forecast Reanalysis]([../Datasets/Datasets.md](https://rda.ucar.edu/datasets/).
For years between 1979 and 2010 the script downloads the data from dataset ds093.1, for more recent years the dataset is ds094.0. The two datasets have different spatial resolutions. 
This create downloads all days beteen Februaray and July and creates one .rds file for each  year. 
`-R`,`-L`,`-S`, and `-N` set the geographic boundaries (see script help for details)

```
for i in {2012..2021}
do 
Rscript windscape_download.R -m 2 -M 7 -d 1 -D 31 -y $i -Y $i -L 345 -R 40 -S 25 -N 60 -o windscape_$i
done
```


The analysis of the data was done as follow. `list_rds_files` is a file with a list of rds files to be combined for the analysis. `samples_coord.csv` is a csv file containing names and coordinates of the samples to be included.  

```
Rscript windscape_analysis.R -L 345 -R 40 -S 25 -N 60 -o windscape_2012-2021 -s samples_coord.csv -l list_rds_files
```
The output consists of two matrixes: 
`windscape_2012-2021.wind_distance_asym.csv` is an asymmetric matrix with each cell (i,j) containing the distance in wind hours moving from individual i towards individual j (from rows to columns). 
`windscape_2012-2021.wind_distance_sym.csv` is an symmetric matrix with each cell containing the average distance in wind hours moving between two individuals in either direction.

The plot of the [wind connectivity surfaces](./) for focal points was produced with the following script:

```
Rscript windscape_plot.R -L 345 -R 40 -S 25 -N 60 -o windscape_2012-2021 -s samples_coord.csv -l list_rds_files -P 0,52,12,53,-5,38,17,41

```

