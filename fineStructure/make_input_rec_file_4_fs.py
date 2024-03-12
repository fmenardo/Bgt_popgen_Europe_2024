import csv
import argparse




parser = argparse.ArgumentParser()

parser.add_argument('-rec','--bp_recombination_map', metavar='', help='input per base recombination map', type =str,nargs=1)
parser.add_argument('-i','--input', metavar='', help='input file with list of positions',type=str,nargs=1)
parser.add_argument('-o','--output', metavar='', help='output stem',type=str,nargs=1)
parser.add_argument('-sep_chr' ,'--separate_chr', default= False, help='separate output files by chromosome', action='store_true')
arguments = parser.parse_args()







# read bp rate file in list
with open(arguments.bp_recombination_map[0], 'r') as csvfile:
	reader = csv.reader(csvfile, delimiter='\t')
	map  = list(reader)



list1=[]
list2=[]

# read list of positions (corresponding o SNPS into list)

with open(arguments.input[0], 'r') as csvfile:
	csv_reader = csv.reader(csvfile)

	count=0
	for row in csv_reader:
		chr, bp = row
		list1.append(row)
		list2.append(bp)




# Append the duplicated column to the list of snp (from the vcf) so that each line corredpond to the interval between snp and following snp
for i, row in enumerate(list1):
	if i==0:
		chr_old= row[0]
	if chr_old != row[0]:
		list1[i-1]=[chr_old,bp_old]
		chr_old=row[0]
	if i <= len(list1)-2:
		row.append(list2[i+1])
	chr_old=row[0]
	bp_old=row[1]


results = {str(i): ["start.pos recom.rate.perbp"] for i in range(1, 12)}


# loop through  list of snp intervals
for i, pairs in enumerate(list1):
	chr=pairs[0]
	n1=pairs[1]
	if len(pairs) > 2:
		n2=pairs[2]
	else:
		n2 = None

	dist=[]
	rec=[]
	dist_tot=int(0)
	rec_tot=float(0)
	flag=0
#	print(pairs)
# loop through list of recobination intervals to calculate average bp recombination rate between the two snps
	for k, row in enumerate(map):
		CHR=row[1]
		rec_r=row[3]
		temp=row[2].split('-')
		start=temp[0]
		end=temp[1]
#		print ("map    "+ str(start)+ "    " + str(end))
		if n2 == None:
			results[str(chr)].append(str(n1)+ " -9")
			break
		if chr == CHR:
			# if the first snp is included in the current iterval
			if int(start) <= int(n1) and int(end) >= int(n1):
				# if also second snp is included in current interval
				if int(n2) <= int(end):
					dist.append(int(n2)-int(n1))
					rec.append(rec_r)
					for t in range(len(dist)):
						dist_tot=dist_tot+dist[t]
						rec_tot=rec_tot+(dist[t]*float(rec[t]))
					final_rec= "{:.10f}".format(rec_tot/dist_tot)
					results[str(chr)].append(str(n1)+ " " + str(final_rec))
					break
				# if second snp is not included
				else:
					dist.append(int(end)-int(n1))
					rec.append(rec_r)
					flag=1
			# if current interval is past first snp, but flag is 1  (therefore we are still recording beacuse not yet arrived at 2nd snp)
			elif int(n1) < int(start) and flag ==1:
					# if second snp is in the current interval
					if int(n2) <= int(end):
						dist.append(int(n2)-int(start))
						rec.append(rec_r)
						for t in range(len(dist)):
							dist_tot=dist_tot+int(dist[t])
							rec_tot=rec_tot+(int(dist[t])*float(rec[t]))
						final_rec= "{:.10f}".format(rec_tot/dist_tot)
						results[str(chr)].append(str(n1)+ " " + str(final_rec))
						break
					# if 2nd snp is still beyon current interval
					elif int(n2) > int(end):
						dist.append(int(end)-int(start))
						rec.append(rec_r)



#print (results)


if(arguments.separate_chr):

	for key, values in results.items():
		with open(arguments.output[0]+'_chr_'+str(key)+'_cp_rec_file.txt', 'w') as file:
			for item in values:
				file.write(str(item) + '\n')

else:
	with open(arguments.output[0]+'_cp_rec_file.txt', 'w') as file:

		for key in results.keys():
			count=0
			for item in results[key]:
				if int(key) == 1:
					file.write(str(item) + '\n')
				elif count >= 1:
					file.write(str(item) + '\n')
				count = count+1
