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

We create a file listing all rds files, and a a csv file containing names and coordinates of the samples to be included (the [Europe+ dataset](../Datasets/Datasets.md)) 

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

Both matrixes have dimensions 401x401, as 14 samples with unknown coordinates have been excluded from the Europe+ dataset (415 samples).

The plot of the [wind connectivity surfaces](./2012-2021_allsites_wind_map.pdf) for some representative focal points was produced with the following script:

```
Rscript windscape_plot.R -L 345 -R 40 -S 25 -N 60 -o windscape_2012-2021 -s samples_coord.csv -l list_wind_files -P 0,52,12,53,-5,38,17,41

```

## Software versions
Many of the spatial packages are in development at the time of writing, and different versions can easily break the code. This analysis was run in these environments:


- Download R sessionInfo()
```
R version 4.2.1 (2022-06-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

locale:
 [1] LC_CTYPE=en_US.UTF-8    LC_NUMERIC=C            LC_TIME=C              
 [4] LC_COLLATE=en_US.UTF-8  LC_MONETARY=C           LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=C              LC_NAME=C               LC_ADDRESS=C           
[10] LC_TELEPHONE=C          LC_MEASUREMENT=C        LC_IDENTIFICATION=C    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ncdf4_1.22      terra_1.7-71    abind_1.4-5     lubridate_1.9.2
 [5] forcats_1.0.0   stringr_1.5.0   dplyr_1.1.2     purrr_1.0.1    
 [9] readr_2.1.4     tidyr_1.3.0     tibble_3.2.1    ggplot2_3.4.2  
[13] tidyverse_2.0.0 raster_3.6-26   sp_2.1-3        windscape_0.1.0
[17] argparser_0.7.1

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.12      pillar_1.9.0     compiler_4.2.1   tools_4.2.1     
 [5] lifecycle_1.0.3  gtable_0.3.3     lattice_0.21-8   timechange_0.2.0
 [9] pkgconfig_2.0.3  rlang_1.1.3      cli_3.6.1        withr_2.5.0     
[13] generics_0.1.3   vctrs_0.6.3      hms_1.1.3        grid_4.2.1      
[17] tidyselect_1.2.0 glue_1.6.2       R6_2.5.1         fansi_1.0.4     
[21] tzdb_0.3.0       magrittr_2.0.3   scales_1.3.0     codetools_0.2-19
[25] colorspace_2.1-0 utf8_1.2.3       stringi_1.7.12   munsell_0.5.0   

```
- Analysis R sessionInfo()

```
R version 4.2.1 (2022-06-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] maps_3.4.1      gdistance_1.6.4 Matrix_1.6-0    igraph_1.5.0   
 [5] lubridate_1.9.2 forcats_1.0.0   stringr_1.5.0   dplyr_1.1.2    
 [9] purrr_1.0.1     readr_2.1.4     tidyr_1.3.0     tibble_3.2.1   
[13] ggplot2_3.4.2   tidyverse_2.0.0 raster_3.6-26   sp_2.1-3       
[17] windscape_0.1.0 argparser_0.7.1

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.12      pillar_1.9.0     compiler_4.2.1   tools_4.2.1     
 [5] lifecycle_1.0.3  gtable_0.3.3     lattice_0.21-8   timechange_0.2.0
 [9] pkgconfig_2.0.3  rlang_1.1.3      cli_3.6.1        terra_1.7-71    
[13] withr_2.5.0      generics_0.1.3   vctrs_0.6.3      hms_1.1.3       
[17] grid_4.2.1       tidyselect_1.2.0 glue_1.6.2       R6_2.5.1        
[21] fansi_1.0.4      tzdb_0.3.0       magrittr_2.0.3   scales_1.3.0    
[25] codetools_0.2-19 geosphere_1.5-14 colorspace_2.1-0 utf8_1.2.3      
[29] stringi_1.7.12   munsell_0.5.0
```
