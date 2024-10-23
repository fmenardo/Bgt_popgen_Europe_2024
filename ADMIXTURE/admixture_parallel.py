## modified Fabrizio's script

import argparse
import os
import subprocess
import random

admixture_path = "/data/admixture"
wd = "/tritici_ALL_25kb_0.1_LDp_admixture/"

def admixture(bed_stem,k,r,n_cores,cv,seed):

	admx= str(admixture_path)+" --cv=" + str(cv) + " -s " + str(seed) + " " + wd + str(bed_stem)+ ".bed --haploid=\"*\" "  + str(k)+ " -j" + str(n_cores) +" | tee logK" + str(k) +".r" + str(r)+ ".out"

	p = subprocess.Popen(admx, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

	OUT, error = p.communicate()
	log.append(admx)
	OUT_log.append(str(OUT))
	error_log.append(str(error))

	os.rename(str(bed_stem)+'.'+str(k)+'.P',str(bed_stem)+'.k'+str(k)+'.r'+str(r)+'.P')
	os.rename(str(bed_stem)+'.'+str(k)+'.Q',str(bed_stem)+'.k'+str(k)+'.r'+str(r)+'.Q')

	return()

# wrapper

parser = argparse.ArgumentParser()
#parser.add_argument('-bp','--bed_path', metavar='', default='', help='path to directory with input bed file', nargs=1,type=str)
parser.add_argument('-B','--bed_stem', metavar='', default='', help='stem of the input bed file', nargs=1,type=str)
parser.add_argument('-K','--number_of_K', metavar='', default=[10], help='number of K (from 1 to K) for the admixture analysis (default=10)', nargs='*',type=int)
parser.add_argument('-r','--number_of_replicates', metavar='',default=[1], help='replicate number of the particular run (for array jobs)', nargs='*',type=int)
parser.add_argument('-c','--number_of_cores', metavar='',default=[1], help='number of cores to use (default=1)', nargs=1,type=int)
parser.add_argument('-o','--out', default='', metavar='',help='stem name for output logs', nargs= 1,type=str)
parser.add_argument('-cv','--n_fold_cross_validation', metavar='', default=[5], help='n-fold cross validation for admixture (default=5)', nargs='*',type=int)

arguments = parser.parse_args()

bed = arguments.bed_stem[0]
K=arguments.number_of_K[0]
R=arguments.number_of_replicates[0]
n_cores=arguments.number_of_cores[0]
out=arguments.out[0]
cv=arguments.n_fold_cross_validation[0]

out_stem=str(out)+"_admixture"

out_path = wd + "r"+str(R)
if not os.path.exists(out_path):
        os.makedirs(out_path)

os.chdir(out_path)

log=[]
error_log=[]
OUT_log=[]

seed = random.randint(0,20000)

# run admixture
for k in range (1, int(K)+1):
	admixture(str(bed),k,R,n_cores,cv,seed)


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


