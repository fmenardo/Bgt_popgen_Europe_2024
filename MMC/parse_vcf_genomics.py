import vcfpy
import argparse



parser = argparse.ArgumentParser()

#parser.add_argument('INFILE',type=str,help='path to the newick tree')
parser.add_argument('-vcf','--input_vcf', metavar='', help='input vcf file' ,default='', type =str,nargs=1)
#parser.add_argument('-K','--number_of_K', metavar='', default=[10], help='number of K (from 1 to K) for the admixture analysis (default=10)', nargs='*',type=int)
#parser.add_argument('-r','--number_of_replicates', metavar='',default=[1], help='how many times to run the admixture analysis for each K (default=1)', nargs='*',type=int)
#parser.add_argument('-anc','--ancestor', metavar='',default='', help='list of outgroups', nargs=1,type=str)
parser.add_argument('-o','--out', default='', metavar='',help='output file', nargs= 1,type=str)
#parser.add_argument('-cv','--n_fold_cross_validation', metavar='', default=[5], help='n-fold cross validation for admixture (default=5)', nargs='*',type=int)



arguments = parser.parse_args()




# Open vcf file, this will read in the header
reader = vcfpy.Reader.from_path(arguments.input_vcf[0])

# Build and print header
header = ['#CHROM', 'POS', 'REF', 'ALT'] + reader.header.samples.names


# initialize fasta file
genome_size=0
OUT={}
Bgt_samples = reader.header.samples.names


for record in Bgt_samples:
	OUT.update({str(record):""})


#print (OUT)


for record in reader:
	alt=[]
	geno = [call.data.get('GT') or './.' for call in record.calls]

	if len(record.ALT) > 0:
		if len(record.ALT[0].value) > 1:
			continue #not a SNP
		if record.ALT[0].value == "*":
			continue # spanning deletion
		for i in range(len(record.ALT)):
			alt.append(record.ALT[i].value)
	if len(record.REF) > 1:
		continue #not a SNP

	if  all(value == "./." for value in geno):  # only missing data
		continue

	ref=record.REF
	genome_size+=1


	for z in range(len(geno)):
		if geno[z] == "0":
			geno[z] = str(ref)
		elif geno[z] == "1":
			geno[z] = str(alt[0])
		elif geno[z] == "2":
			geno[z] = str(alt[1])
		elif geno[z] == "3":
			geno[z] = str(alt[2])
		elif geno[z] == "./.":
			geno[z] = "-"
	pairs = zip(reader.header.samples.names, geno)

	my_dict = dict(pairs)

	for name in Bgt_samples:
		OUT[name] += my_dict[name]


print(genome_size)
f = open(arguments.out[0], "w")
for name_out, seq_out in OUT.items():
		f.write(">" + name_out + "\n")
		f.write(seq_out+"\n")
f.close()

f = open(str(arguments.out[0])+".genome_size", "w")

f.write(str(genome_size))

f.close()

