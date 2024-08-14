import vcfpy
import argparse



parser = argparse.ArgumentParser()

#parser.add_argument('INFILE',type=str,help='path to the newick tree')
parser.add_argument('-vcf','--input_vcf', metavar='', help='input vcf file' ,default='', type =str,nargs=1)
#parser.add_argument('-K','--number_of_K', metavar='', default=[10], help='number of K (from 1 to K) for the admixture analysis (default=10)', nargs='*',type=int)
#parser.add_argument('-r','--number_of_replicates', metavar='',default=[1], help='how many times to run the admixture analysis for each K (default=1)', nargs='*',type=int)
parser.add_argument('-anc','--ancestor', metavar='',default='', help='list of outgroups', nargs=1,type=str)
parser.add_argument('-o','--out', default='', metavar='',help='output file', nargs= 1,type=str)
#parser.add_argument('-cv','--n_fold_cross_validation', metavar='', default=[5], help='n-fold cross validation for admixture (default=5)', nargs='*',type=int)



arguments = parser.parse_args()



## open list of outgroups and put names i list
with open(arguments.ancestor[0], 'r') as file:
    Bgs_samples = [line.rstrip() for line in file]

#print(Bgs_samples)



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



# initialize fasta file
genome_size=0
OUT={}
Bgt_samples = [reader.header.samples.names[i] for i in Bgt_ind]


OUT.update({"ANC":""})

for record in Bgt_samples:
	OUT.update({str(record):""})


#print (OUT)


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
		if record.ALT[0].value == "*":
			print (record.ALT)
			continue # spanning deletion

	if len(record.REF) > 1:
		continue #not a SNP
	count_Bgs=0
	if "1" in Bgs:
		count_Bgs+=2
	if "0" in Bgs:
		count_Bgs+=1


	if count_Bgs == 0:
		continue #Bgs has only missing data
	if count_Bgs == 3:
		continue #Bgs polymorphic
	else:            #Bgs monomorphic

		count_Bgt=0

		if "1" in Bgt:
			count_Bgt+=2
		if "0" in Bgt:
			count_Bgt+=1

		if count_Bgt == 2:
			genome_size+=1  # Bgt monomorphic
			continue

		if count_Bgt == 1:
			genome_size+=1   # Bgt_monomorphic
			continue


		if count_Bgt == 3:      # Bgt polymorphic
			genome_size+=1
			ref=record.REF
			alt=record.ALT[0].value    # works only for biallelic snps

			for z in range(len(geno)):
				if geno[z] == "0":
					geno[z] = str(ref)
				elif geno[z] == "1":
					geno[z] = str(alt)


			pairs = zip(reader.header.samples.names, geno)

			my_dict = dict(pairs)

			for name in Bgt_samples:
				OUT[name] += my_dict[name]

			if count_Bgs ==1:
				OUT["ANC"] += str(ref)
			elif count_Bgs ==2:
				OUT["ANC"] += str(alt)




print(genome_size)
f = open(arguments.out[0], "w")
for name_out, seq_out in OUT.items():
		f.write(">" + name_out + "\n")
		f.write(seq_out+"\n")
f.close()

f = open(str(arguments.out[0])+".genome_size", "w")

f.write(str(genome_size))

f.close()

