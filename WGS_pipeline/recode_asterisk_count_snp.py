"""
Removes sites with "*" in ALT column (spanning deletions) from a snp vcf file.
"""

import argparse
import subprocess
import os
import re
import gzip

os.chdir("../data")

def recode_vcf(in_file,out_file):
    bunchsize = 1000000
    bunch = []
    ast = 0
    snp = 0
    
    with gzip.open(f"{str(in_file)}", mode = 'rt') as file, open(f"{str(out_file)}.vcf", mode = 'a+') as F:
        for line in file:
            line=line.strip()
            if line.startswith("#"):        # header
                bunch.append(line)
                
            else:
                chunks = re.split('\s+',line)
                alt = re.split(',',chunks[4])
                if "*" in alt:
                    ast += 1
                    if len(alt) == 2:
                        if len(alt[0]) == 1 and len(alt[1]) == 1:
                            snp += 1  
                else:
                    bunch.append(line)
                

            if len(bunch) == bunchsize:
                F.write('\n'.join(bunch))
                F.write('\n')
                bunch = []
        F.write('\n'.join(bunch))
        F.write('\n')
    return ast, snp    
        
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', metavar = '', help='name of compressed input vcf file to be recoded')
parser.add_argument('-o', '--output', metavar='',help='prefix for output vcf file' )

ag = parser.parse_args()

if __name__== '__main__':
    myfunc = recode_vcf(ag.input,ag.output)
    print(f"removed {myfunc[0]} sites of which {myfunc[1]} were sites with only biallelic SNPs.")
    
