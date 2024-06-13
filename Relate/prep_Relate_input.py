import vcfpy
import argparse
import csv


parser = argparse.ArgumentParser()

#parser.add_argument('INFILE',type=str,help='path to the newick tree')
parser.add_argument('-vcf','--input_vcf', metavar='', help='input vcf file' ,default='', type =str,nargs=1)
#parser.add_argument('-K','--number_of_K', metavar='', default=[10], help='number of K (from 1 to K) for the admixture analysis (default=10)', nargs='*',type=int)
parser.add_argument('-meta','--metainfo', metavar='', help='path to csv file with metainformation', nargs=1,type=str)
parser.add_argument('-anc','--ancestor', metavar='',default='', help='list of outgroups', nargs=1,type=str)
parser.add_argument('-o','--out', default='', metavar='',help='output file', nargs= 1,type=str)
parser.add_argument('-onc','--out_not_chr', default='', metavar='',help='output file not specific for one chromosome (.sample, .poplabel)', nargs= 1,type=str)

#parser.add_argument('-cv','--n_fold_cross_validation', metavar='', default=[5], help='n-fold cross validation for admixture (default=5)', nargs='*',type=int)



arguments = parser.parse_args()



## open list of outgroups and put names i list
with open(arguments.ancestor[0], 'r') as file:
    Bgs_samples = [line.rstrip() for line in file]


# Open vcf file, this will read in the header
reader = vcfpy.Reader.from_path(arguments.input_vcf[0])

# Build and print header
header = ['#CHROM', 'POS', 'REF', 'ALT'] + reader.header.samples.names



#find indeces of outgoups and of ingroup

Bgs_ind=[]
Bgt_ind=[]
for i in range(len(reader.header.samples.names)):
#	print (i)
	flag=0
	for sample in Bgs_samples:
		if sample == reader.header.samples.names[i]:
			flag=1
	if flag==1:
		Bgs_ind.append(i)
	else:
		Bgt_ind.append(i)



# make .sample file
Bgt_samples = [reader.header.samples.names[i] for i in Bgt_ind]
sample_file =[]
sample_file.append("ID_1 ID_2 missing")
sample_file.append("0 0 0")


for record in Bgt_samples:
	sample_file.append(str(record)+" NA 0")


f = open(str(arguments.out_not_chr[0])+".sample", "w")


for item in sample_file:
	f.write(str(item) + "\n")

f.close()

# make .poplabels file
rowcount=0
sample_count=0
poplab_file =[]
poplab_file.append("sample population group sex")
with open(arguments.metainfo[0], newline='') as csvfile:

	csv_reader = csv.reader(csvfile)

	selected_col = ['Sample.Name', 'fs_level_4', 'fs_level_10']

	for values  in csv_reader:
		rowcount += 1
		if rowcount ==1:
			indexes = [values.index(col) for col in selected_col]
		elif rowcount > 1:
			if values[indexes[0]] in Bgt_samples:
				sample_count +=1
				poplab_file.append(str(values[indexes[0]]) + " " + str(values[indexes[1]]) + " " + str(values[indexes[2]])+ " 1")


f = open(str(arguments.out_not_chr[0])+".poplabels", "w")


for item in poplab_file:
        f.write(str(item) + "\n")

f.close()



### main loop
dist=0
count_SNP=0
hap_file=[]
dist_file=[]
dist_file.append("#pos dist")
genome_size=0
for record in reader:
	geno = [call.data.get('GT') or './.' for call in record.calls]
	list_geno = list(geno)
	Bgs = [geno[i] for i in Bgs_ind]
	Bgt = [geno[i] for i in Bgt_ind]

	if "./." in Bgt:
		continue #Bgt has at least 1 missing data
	if len(record.ALT) > 0:
		if len(record.ALT[0].value) > 1:
			continue #not a SNP
	if len(record.REF) > 1:
		continue #not a SNP
	count_Bgs=0
	if "1" in Bgs:
		count_Bgs+=2
	if "0" in Bgs:
		count_Bgs+=1


	if count_Bgs == 0:
		continue #Bgs has only missing data
	elif count_Bgs == 3:
		continue #Bgs polymorphic
	else:            #Bgs monomorphic
		dist+=1
		genome_size=0
		count_Bgt=0

		if "1" in Bgt:
			count_Bgt+=2
		if "0" in Bgt:
			count_Bgt+=1

		if count_Bgt == 2: # Bgt monomorphic
			continue

		elif count_Bgt == 1: # Bgt_monomorphic
			continue

		elif count_Bgt == 3:      # Bgt polymorphic
			count_SNP+=1
			if count_SNP == 1:
				dist=0
			elif count_SNP > 1:
				dist_file.append(str(pos_old) + " " + str(dist))
				dist=0
			ref=record.REF
			alt=record.ALT[0].value    # works only for biallelic snps
							# ancestral is 0
			if count_Bgs ==1:
				ANC = str(ref)
				DER = str(alt)
							# ancestral is 1 => invert all genotypes
			if count_Bgs ==2:
				ANC = str(alt)
				DER = str(ref)
				for i in range(len(Bgt)):
					if Bgt[i] == 0:
						Bgt[i] = 1
					elif Bgt[i] == 1:
						Bgt[i] = 0
			str_geno = ' '.join(map(str, Bgt))

			line = (str(record.CHROM)+ " SNP" + str(count_SNP)+ " " + str(record.POS) + " " + ANC + " " + DER + " "  +str_geno)
			hap_file.append(line)
			pos_old=record.POS

dist_file.append(str(pos_old) + " " + str(dist))

# write hap file


f = open(str(arguments.out[0])+".haps", "w")


for item in hap_file:
        f.write(str(item) + "\n")
f.close()


# write dist file


f = open(str(arguments.out[0])+".dist", "w")


for item in dist_file:
        f.write(str(item) + "\n")
f.close()

print(genome_size)
