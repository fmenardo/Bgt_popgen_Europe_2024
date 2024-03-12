import csv
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('-rec','--recombination_map', metavar='', help='input recombination map', type =str,nargs=1)
parser.add_argument('-o','--output', metavar='', help='output stem',type=str,nargs=1)
parser.add_argument('-uni' ,'--uniform', default= False, help='make a uniform recombination map, for each chromodome take the total cM and divide by the physical length (also create a unifomr map over the whole genome)', action='store_true')


arguments = parser.parse_args()



results = []
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

chromosome_cM = {str(i): [None] for i in range(1, 12)}


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
			results.pop()
			rec_r= (float(cM_old)-float(cM_start_old))/(int(chromosome_lengths[str(chr_old)])-int(bp_start_old))/100
			results.append("chr\t" + str(chr_old)+ "\t" +str(bp_start_old)+"-"+str(chromosome_lengths[str(chr_old)]) + "\t" + str("{:.10f}".format(rec_r)))
			chr_old=chr
			cM_old=cM
			bp_old=bp
			continue


		# if there is the cM count changes along the map (new recombination)
		if float(cM_start) != float(cM):
			if new_chr == 1:
				rec_r= float(cM)/(int(bp))/100
				new_chr=2
				results.append("chr\t" + str(chr) + "\t"+ "1 -"+str(bp) + "\t" + str("{:.10f}".format(rec_r)))
			else:
				rec_r= (float(cM)-float(cM_start))/(int(bp)-int(bp_start))/100
				results.append("chr\t" + str(chr) +"\t"+ str(bp_start)+"-"+str(bp) + "\t" + str("{:.10f}".format(rec_r)))

			cM_start_old=cM_start
			cM_start=cM

			bp_start_old=bp_start
			bp_start=int(bp)

		cM_old=cM
		bp_old=bp
		chr_old=chr


results.pop()

rec_r= (float(cM)-float(cM_start_old))/(int(chromosome_lengths[chr])-int(bp_start_old))/100

results.append("chr\t" + str(chr) +"\t"+str(bp_start_old)+"-"+str(chromosome_lengths[chr]) + "\t" + str("{:.10f}".format(rec_r)))

chromosome_cM[chr]=cM_old



with open(str(arguments.output[0])+'_bp_recombination_rates.txt', 'w') as file:
	for item in results:
		file.write(str(item) + '\n')


if arguments.uniform:
	results=[]
	for i in  chromosome_lengths:
		results.append("chr\t" + str(i) +"\t1-"+str(chromosome_lengths[i]) +"\t" + str("{:.10f}".format((float(chromosome_cM[i])/int(chromosome_lengths[i]))/100)))


	with open(str(arguments.output[0])+'_bp_recombination_rates_uniform_o_chr.txt', 'w') as file:
		for item in results:
			file.write(str(item) + '\n')

	results=[]
	tot_cM=0
	genome_l=0
	for i in  chromosome_lengths:
		tot_cM= tot_cM+float(chromosome_cM[i])
		genome_l= genome_l+int(chromosome_lengths[i])
	for i in  chromosome_lengths:
		results.append("chr\t" + str(i) +"\t1-"+str(chromosome_lengths[i]) +"\t" + str("{:.10f}".format((tot_cM/genome_l)/100)))


	with open(str(arguments.output[0])+'_bp_recombination_rates_uniform_o_genome.txt', 'w') as file:
		for item in results:
			file.write(str(item) + '\n')


