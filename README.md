# Population genomics and molecular epidemiology of wheat powdery mildew in Europe

Here you find the code to reproduce the analyses performed in **Population genomics and molecular epidemiology of wheat powdery mildew in Europe** (Jigisha et al. 2024)

With this repository we try to ensure reproducibility of all the analyses performed for this study. While we strived to document our code and facilitate the reproduction of our analyses, this is not a ready to run repository. At the very least you will need to install the correct software and adjust the paths to be able to run this code.

The reference genome assembly as well as the VCF file used in this study  are available [here]{https://doi.org/10.5281/zenodo.13903934}.

- Data preparation
  - [WGS pipeline](WGS_pipeline/WGS_pipeline.md)
  - [Datasets](Datasets/Datasets.md)
  - [Recombination map](recombination_map/recombination_map.md)
  - [Distance matrix](distance_matrix/distance_matrix.md)
  - [Windscape](windscape/windscape.md)
- Population structure and spatial population genetics
  - [Principal component analysis](PCA/PCA.md) 
  - [ADMIXTURE](ADMIXTURE/ADMIXTURE.md)
  - [fineStructure](fineStructure/fineStructure.md)
  - [Isolation by Distance](Isolation_by_distance/Isolation_by_distance.md)
  - [FEEMS](FEEMS/FEEMS.md)
- Demographic inference
   - [Demographic inference](MMC/MMC.md)
- Genome scans and summary statistics
  - [Summary statistics](summary_statistics/summary_statistics.md)
  - [isoRelate](isoRelate/isoRelate.md)
  - [LD decay](Linkage_Disequilibrium/LD.md)
- AvrPm17
  - [Haplotypes](AvrPm17_haplotypes/AvrPm17_haplotypes.md) 
  - [Identity by descent clusters](AvrPm17_isoRelate/AvrPm17_isoRelate.md)
  - [Quantification of HR intensity](AvrPm17_infiltration/AvrPm17_infiltration.md)
- Redundancy analysis
  - [RDA](RDA/RDA.md) 
