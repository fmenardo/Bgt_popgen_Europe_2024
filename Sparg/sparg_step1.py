import argparse
import tskit
import subprocess



parser = argparse.ArgumentParser()

parser.add_argument('-in','--input', metavar='', help='input stem without ".anc" and ".mut"', type =str,nargs=1)
parser.add_argument('-o','--output', metavar='', help='path to output folder + stem',type=str,nargs=1)
parser.add_argument('-coal','--coal', metavar='', help='.coal file', type =str,nargs=1)
parser.add_argument('-chr','--chr', metavar='', help='chromosome number ',type=int,nargs=1)
parser.add_argument('-mut','--mut', metavar='', help='mutation rate ',type=str,nargs=1)
parser.add_argument('-ns','--ns', metavar='', help='number of times branches are resampled by relate ',type=int,nargs=1)
parser.add_argument('-min_snps','--min_snps', metavar='', help='min number of snps in a tree to be retained ',type=int,nargs=1)
parser.add_argument('-min_win','--min_win', metavar='', help='min distance between trees (in bp)  ',type=int,nargs=1)
parser.add_argument('-nth_tree' ,'--nth_tree', default= False, help='pick every 250th tree instead of the trees with more than X mutations', action='store_true')


arguments = parser.parse_args()



trees=[]
ixs_start=[]
ixs_end=[]
min_mut= arguments.min_snps[0]
min_window= arguments.min_win[0]
nsamples = str(arguments.ns[0])
m= arguments.mut [0]
chr = str(arguments.chr[0])
output= arguments.output[0]                                                                  

coal= arguments.coal[0]

input= arguments.input[0]
dist= input+".dist"



inf_ts = tskit.load(input+".trees")


if arguments.nth_tree:

	n=len(inf_ts.trees())

	trees = [i for i in range(0,n+1,250)] #which trees to sample

	for tree in inf_ts.trees():
		if tree.index in trees:
			ixs_start.append(int(tree.interval.left))
			ixs_end.append(int(tree.interval.right))

else:


	end_old= -(2*min_window)  # trick for the first record
	for tree in inf_ts.trees():
		mutation_count = sum(1 for site in tree.sites() for _ in site.mutations)
		if mutation_count >= min_mut:
			if ((int(tree.interval.left) - end_old) > min_window):
				trees.append(tree.index)
				ixs_start.append(tree.interval.left)
				ixs_end.append(tree.interval.right)
				end_old = int(tree.interval.right)

# run relate to resample branch lengths "nsamples" times
for i in range(len(trees)):

	command = "~/data/relate/scripts/SampleBranchLengths/SampleBranchLengths.sh " + \
		" -i " + input + \
		" --dist "+ dist + \
		" --coal "+ coal + \
		" -o " + output + ".chr" + chr+ "." + str(i) +\
		" -m " + m + \
		" --format n --num_samples " + str(nsamples) + \
		" --first_bp "+ str(int(ixs_start[i])) + \
		" --last_bp " + str(int(ixs_end[i])) + \
		" --seed 564564"



	p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	out, error = p.communicate()


#import sys
#sys.path.append("~/projects/project_tritici_fabrizio/sparg/")
#from newick import from_newick
#from utils import get_shared_times
#import numpy as np
