# Windscape
First the wind data was downloaded from the [Climate System Forecast Reanalysis]([../Datasets/Datasets.md](https://rda.ucar.edu/datasets/).
For years between 1979 and 2010 the script downloads the data from dataset ds093.1, for more recent years the dataset is ds094.0. The two datasets have different spatial resolutions. 
This create downloads all days beteen Februaray and July and creates one .rds file for each  year.

```
for i in {2012..2021}
do 
Rscript windscape_download.R -m 2 -M 7 -d 1 -D 31 -y $i -Y $i -L 345 -R 40 -S 25 -N 60 -o windscape_2011
done
```


The analysis of the data was done withthe following script, where 


```
Rscript windscape_analysis.R -L 345 -R 40 -S 25 -N 60 -o windscape_2012-2021 -s samples_coord.csv -l list_rds_files
```

