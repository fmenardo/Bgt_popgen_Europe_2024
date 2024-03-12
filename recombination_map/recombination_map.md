

# extract info from genetic maps 
head -n 1 ../data/GeneticMap_final_13052018.csv | sed 's/Isolate,//g' | sed 's/,/\n/g' > list_snp_ig_G
head -n 2 ../data/GeneticMap_final_13052018.csv| tail -n 1 | sed 's/^,1/1/g' | sed 's/,/\n/g'  > list_chr_G
head -n 3 ../data/GeneticMap_final_13052018.csv| tail -n 1 | sed 's/^,0/0/g' | sed 's/,/\n/g'  > list_rec_G

#extract physical positions to attach to genetic map
grep -w -f list_snp_ig_G ../data/Positions_THUN12X96224.csv  > Positions_THUN12X96224_present_in_G.csv

cut -d ";" -f 3 Positions_THUN12X96224_present_in_G.csv  > list_physical_position_G

paste  -d "," list_snp_ig_G list_chr_G list_rec_G list_physical_position_G  > THUN12x96224_genetic_map_in_cM_+_phy_distance
