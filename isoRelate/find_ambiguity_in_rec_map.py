import csv
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('-rec','--recombination_map', metavar='', help='input recombination map with physical postions anad cM format: 1,2,3,4 wher 1=snp id , 2= chromosome number, 3 = cM, 4 = bp', type =str,nargs=1)


arguments = parser.parse_args()



## thin recombination map

thin_rec_map=[]
intervals=[]
with open(arguments.recombination_map[0], 'r') as csvfile:
	csv_reader = csv.reader(csvfile)
	row_old=""
	CHR_OLD=1
	CM_OLD=0
	count=0
	BP_OLD=0
	for row in csv_reader:
		count=count+1
		snpID, CHR, CM, BP = row
		if int(CHR) != int(CHR_OLD):
			thin_rec_map.append(row_old)
			thin_rec_map.append(row)
			CHR_OLD=CHR
			CM_OLD=CM
			row_old=row
			BP_OLD=BP
			continue
		elif ((int(CHR) == int(CHR_OLD)) and (float(CM) > float(CM_OLD))):
			thin_rec_map.append(row_old)
			thin_rec_map.append(row)
			CHR_OLD=CHR
			CM_OLD=CM
			BP_OLD=BP
		row_old=row

with open(arguments.recombination_map[0]+'.thinned', 'w') as file:
    for item in thin_rec_map:
        file.write(str(item) + '\n')


CHR_OLD=1
BP_OLD=0
CM_OLD=0


list_ctg=["LR026984.1_chr1","LR026985.1_chr2","LR026986.1_chr3","LR026987.1_chr4","LR026988.1_chr5","LR026989.1_chr6","LR026990.1_chr7","LR026991.1_chr8","LR026992.1_chr9","LR026993.1_chr10","LR026994.1_chr11"]

for row in thin_rec_map:
	count=count+1
	snpID, CHR, CM, BP = row
	if ((int(CHR) == int(CHR_OLD)) and (float(CM) > float(CM_OLD))):
		intervals.append(list_ctg[int(CHR)-1]+":"+str(int(BP_OLD)+1) +"-"+str(int(BP)-1))
	CHR_OLD=CHR
	CM_OLD=CM
	BP_OLD=BP


with open(arguments.recombination_map[0]+'.ambiguous_intervals.list', 'w') as file:
    for item in intervals:
        file.write(str(item) + '\n')
