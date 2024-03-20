## modified Fabrizio's script

import argparse
import os
#import csv
#import sys
#import re
import subprocess
#import random
#import shutil

plink_path = "/home/jjigis/data/plink"
admixture_path = "/home/jjigis/data/dist/admixture_linux-1.3.0/admixture"

def make_bed(vcf_file,n_cores):


	zcat= "zcat " +  str(vcf)+ " | sed s/LR0269[0-9][0-9].1_chr//g > "+ str(out) + "_renamed_chr.vcf"

	p = subprocess.Popen(zcat, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

	OUT, error = p.communicate()
	log.append(zcat)
	OUT_log.append(str(OUT))
	error_log.append(str(error))

	#gzip= "gzip -f " +  str(out) + "_renamed_chr.vcf"
	gzip= "bgzip " +  str(out) + "_renamed_chr.vcf"

	p = subprocess.Popen(gzip, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

	OUT, error = p.communicate()
	log.append(gzip)
	OUT_log.append(str(OUT))
	error_log.append(str(error))


	plink= str(plink_path)+" --vcf " +  str(out) + "_renamed_chr.vcf.gz --make-bed --allow-extra-chr --double-id --threads " + str(n_cores)+ " --out " + str(out)

	p = subprocess.Popen(plink, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

	OUT, error = p.communicate()
	log.append(plink)
	OUT_log.append(str(OUT))
	error_log.append(str(error))
	return()


def admixture(bed_file,k,r,n_cores,cv):

	admx= str(admixture_path)+" --cv=" +str(cv) + " " + str(bed_file)+ " --haploid=\"*\" "  + str(k)+ " -j" + str(n_cores) +" | tee logK" + str(k) +".r" + str(r)+ ".out"

	p = subprocess.Popen(admx, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

	OUT, error = p.communicate()
	log.append(admx)
	OUT_log.append(str(OUT))
	error_log.append(str(error))

	os.rename(str(out)+'.'+str(k)+'.P',str(out)+'.k'+str(k)+'.r'+str(r)+'.P')
	os.rename(str(out)+'.'+str(k)+'.Q',str(out)+'.k'+str(k)+'.r'+str(r)+'.Q')

	return()

# wrapper

parser = argparse.ArgumentParser()

#parser.add_argument('INFILE',type=str,help='path to the newick tree')
parser.add_argument('-vcf','--input_vcf', metavar='', help='input vcf file' ,default='', type =str,nargs=1)
parser.add_argument('-K','--number_of_K', metavar='', default=[10], help='number of K (from 1 to K) for the admixture analysis (default=10)', nargs='*',type=int)
parser.add_argument('-r','--number_of_replicates', metavar='',default=[1], help='how many times to run the admixture analysis for each K (default=1)', nargs='*',type=int)
parser.add_argument('-c','--number_of_cores', metavar='',default=[1], help='number of cores to use (default=1)', nargs=1,type=int)
parser.add_argument('-o','--out', default='', metavar='',help='stem name for output', nargs= 1,type=str)
parser.add_argument('-cv','--n_fold_cross_validation', metavar='', default=[5], help='n-fold cross validation for admixture (default=5)', nargs='*',type=int)



arguments = parser.parse_args()

vcf=arguments.input_vcf[0]

'''
if vcf[0] == "/":
	path_fixer=''
else:
	#path_fixer='../'
	path_fixer = '../../../data/'
'''
n_cores=arguments.number_of_cores[0]

K=arguments.number_of_K[0]
out=arguments.out[0]
R=arguments.number_of_replicates[0]
cv=arguments.n_fold_cross_validation[0]

out_stem=str(out)+"_admixture"


# create folder 
if not os.path.exists(out_stem):
	os.makedirs(out_stem)
os.chdir(out_stem)
os.system("rm *.out")
os.system("rm *.Q")
os.system("rm *.P")
os.system("rm *.log")




log=[]
error_log=[]
OUT_log=[]

# convert vcf to bed, changing chromosome names

make_bed(vcf,n_cores)


# run admixture
for r in range (1, int(R)+1):
	for k in range (1, int(K)+1):
		admixture(str(out)+".bed",k,r,n_cores,cv)


# make logfiles

logfile= str(out_stem)+".log"

with open(logfile, 'w') as fp:
	fp.write('\n'.join(log))


if len(error_log) > 0:

	logfile= str(out_stem)+".error.log"

	with open(logfile, 'w') as fp:
		fp.write('\n'.join(error_log))



logfile= str(out_stem)+".out.log"

with open(logfile, 'w') as fp:
	fp.write('\n'.join(OUT_log))


grep= "grep -h CV log*.out > " + str(out) + "_admixture_K.log"

p = subprocess.Popen(grep, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

OUT, error = p.communicate()
log.append(grep)
OUT_log.append(str(OUT))
error_log.append(str(error))
