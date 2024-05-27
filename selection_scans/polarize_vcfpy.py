import vcfpy
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-vcf', '--input_vcf', metavar='', help='input vcf file', default='', type=str, nargs=1)
parser.add_argument('-anc', '--ancestor', metavar='', default='', help='list of outgroups', nargs=1, type=str)
parser.add_argument('-o', '--out', default='', metavar='', help='output file', nargs=1, type=str)

arguments = parser.parse_args()

# Open list of outgroups and put names in a list
with open(arguments.ancestor[0], 'r') as file:
    Bgs_samples = [line.rstrip() for line in file]

# Open VCF file
reader = vcfpy.Reader.from_path(arguments.input_vcf[0])

# Modify header to include AA INFO field
reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'AA'), ('Number', '1'), ('Type', 'String'), ('Description', 'Ancestral Allele')]))

# Open output VCF file for writing
writer = vcfpy.Writer.from_path(arguments.out[0], reader.header)

# Find indices of outgroups and ingroup
Bgs_ind = []
Bgt_ind = []
for i, sample_name in enumerate(reader.header.samples.names):
    if sample_name in Bgs_samples:
        Bgs_ind.append(i)
    else:
        Bgt_ind.append(i)

for record in reader:
     line = [record.CHROM, record.POS, record.REF]
     line += [alt.value for alt in record.ALT]
     geno = [call.data.get('GT') or './.' for call in record.calls]
     list_geno = list(geno)
     Bgs = [geno[i] for i in Bgs_ind]
     Bgt = [geno[i] for i in Bgt_ind]
     count_Bgs=0
     if "1" in Bgs:
         count_Bgs+=2
     if "0" in Bgs:
         count_Bgs+=1
     ref=record.REF
     alt=record.ALT[0].value # works only for biallelic snps
     ancestral_allele = ''
     # add "AA=" filter to INFO field
     if count_Bgs == 1:
         ancestral_allele = ref
     elif count_Bgs == 2:
         ancestral_allele = alt
     elif count_Bgs == 3:
         ancestral_allele = f"{record.REF}{record.ALT[0].value}"  # If both REF and ALT alleles are ancestral
     # Add AA information to the INFO field
     if ancestral_allele:
         record.INFO['AA'] = ancestral_allele
     writer.write_record(record)
writer.close()

# Process each record
#for record in reader:
#    # Calculate ancestral allele information
#    ancestral_allele = ''
#    count_Bgs = sum(1 for i in Bgs_ind if record.calls[i].data.get('GT') in ['0/0', '0|0', '0'])
#    if count_Bgs == 1:
#        ancestral_allele = record.REF
#    elif count_Bgs == 2:
#        ancestral_allele = record.ALT[0].value
#    elif count_Bgs == 3:
#        ancestral_allele = f"{record.REF}{record.ALT[0].value}"  # If both REF and ALT alleles are ancestral
#    
#    # Add AA information to the INFO field
#    if ancestral_allele:
#        record.INFO['AA'] = ancestral_allele
#
#    # Write the modified record to the output VCF file
#    writer.write_record(record)
#
## Close the writer
#writer.close()
