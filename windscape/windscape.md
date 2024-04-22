# Windscape
We used the package [windscape](https://github.com/matthewkling/windscape) to infer wind connectivity between sampling locations.

The wind data (u and v comoponent at 10 m above sea level) was downloaded from the [Climate System Forecast Reanalysis](https://rda.ucar.edu/datasets/).
For years between 1979 and 2010 the script downloads the data from dataset ds093.1, for more recent years the dataset is ds094.0. The two datasets have different spatial resolutions. 
The following command downloads the data for all days between February and July and creates one .rds file for each year between 2012 and 2021. See script help (`-h`) for details about arguments. 

```
for i in {2012..2021}
do 
Rscript windscape_download.R -m 2 -M 7 -d 1 -D 31 -y $i -Y $i -L 345 -R 40 -S 25 -N 60 -o windscape_$i
done
```

We create a file listing all rds files, and a a csv file containing names and coordinates of the samples to be included 

```
ls | grep "rds" > list_wind_files

head -n 1 ../Datasets/2022+before2022+2023+ncsu_metadata+fs+admxK7_19032024.csv | cut -f 1,8,9 -d","  > samples_coord.csv
cut -f 1,8,9 -d"," ../Datasets/2022+before2022+2023+ncsu_metadata+fs+admxK7_19032024.csv | grep -f ../Datasets/tritici_extended_europe_2022+before2022+2023+ncsu.args >> samples_coord.csv

```

We calculate distances: 

```
Rscript windscape_analysis.R -L 345 -R 40 -S 25 -N 60 -o windscape_2012-2021 -s samples_coord.csv -l list_wind_files
```

The output consists of two matrixes: 
`windscape_2012-2021.wind_distance_asym.csv` is an asymmetric matrix with each cell (i,j) containing the distance in wind hours moving from individual i towards individual j (from rows to columns). 
`windscape_2012-2021.wind_distance_sym.csv` is an symmetric matrix with each cell containing the average distance in wind hours moving between two individuals in either direction.

The plot of the [wind connectivity surfaces](./2012-2021_allsites_wind_map.pdf) for some representative focal points was produced with the following script:

```
Rscript windscape_plot.R -L 345 -R 40 -S 25 -N 60 -o windscape_2012-2021 -s samples_coord.csv -l list_wind_files -P 0,52,12,53,-5,38,17,41

```

