import csv
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('-rec','--recombination_map', metavar='', help='input recombination map', type =str,nargs=1)
parser.add_argument('-o','--output', metavar='', help='output stem',type=str,nargs=1)
#parser.add_argument('-uni' ,'--uniform', default= False, help='make a uniform recombination map, for each chromodome take the total cM and divide by the physical length (also create a unifomr map over the whole genome)', action='store_true')


arguments = parser.parse_args()



bp_start=1
chr_old=1
cM_start=0
new_chr=1
chromosome_lengths = {
        '1': 12010707,
        '2': 15587074,
        '3': 5931322,
        '4': 11072032,
        '5': 19486627,
        '6': 11625509,
        '7': 16229997,
        '8': 11902262,
        '9': 19721624,
        '10': 10259397,
        '11': 3575586
}


chromosome_names = {
        '1': "LR026984.1_chr1",
        '2': "LR026985.1_chr2",
        '3': "LR026986.1_chr3",
        '4': "LR026987.1_chr4",
        '5': "LR026988.1_chr5",
        '6': "LR026989.1_chr6",
        '7': "LR026990.1_chr7",
        '8': "LR026991.1_chr8",
        '9': "LR026992.1_chr9",
        '10': "LR026993.1_chr10",
        '11': "LR026994.1_chr11"
}



chromosome_cM = {str(i): [None] for i in range(1, 12)}


results = {}
for i in range(1, 12):
    results[i] = ["pos COMBINED_rate Genetic_Map"]


# loop through the recombination map
with open(arguments.recombination_map[0], 'r') as csvfile:
	csv_reader = csv.reader(csvfile)

	for row in csv_reader:
		snpID, chr, cM, bp = row
		# when the chromosome change  I rewrite the last line to extend the interval up to the end of teh chromosome
		if int(chr) != int(chr_old):
			chromosome_cM[chr_old] = cM_old
			cM_start=0
			new_chr = 1
			rec_r= (float(cM_old)-float(cM_start_old))/(int(chromosome_lengths[str(chr_old)])-int(bp_start_old))*1000000
			results[int(chr_old)].append(str(bp_start)+ " " + str(rec_r)+ " " + str(cM_old))
			chr_old=chr
			cM_old=cM
			bp_old=bp
			continue


		# if there is the cM count changes along the map (new recombination)
		if float(cM_start) != float(cM):
			if new_chr == 1:
				rec_r= float(cM)/(int(bp))*1000000
				new_chr=2
				results[int(chr)].append("1 " + str(rec_r) + " 0")
			else:
				rec_r= (float(cM)-float(cM_start))/(int(bp)-int(bp_start))*1000000
				results[int(chr)].append(str(bp_start) + " " + str(rec_r) + " " + str(cM_start))

			cM_start_old=cM_start
			cM_start=cM

			bp_start_old=bp_start
			bp_start=int(bp)

		cM_old=cM
		bp_old=bp
		chr_old=chr



rec_r= (float(cM_old)-float(cM_start_old))/(int(chromosome_lengths[str(chr_old)])-int(bp_start_old))*1000000
results[int(chr)].append(str(bp_start)+ " " + str(rec_r)+ " " + str(cM_old))

chromosome_cM[chr]=cM_old



for key, value in results.items():
    file_name = str(arguments.output[0])+ "_chr" + str(key) +'_Mb_recombination_rates.txt'

    with open(file_name, "w") as file:
        for item in value:
            file.write(f"{item}\n")





