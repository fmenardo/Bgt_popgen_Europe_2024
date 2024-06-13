from tsnewick import from_newick
import numpy as np
import argparse
import tskit
from utils import get_shared_times, chop_shared_times, center_shared_times
import glob
import os
from spacetrees import _log_coal_density
from joblib import Parallel, delayed
import time

# main function
def process_newick(file,tCutoff,coal):
	stss = []
	ctss = []
	output=file.replace(".newick","")
	count=0
	with open(file, mode='r') as f:
		next(f) #skip header
		for line in f: #for each tree sampled at a locus
			count+=1
 			# import
			string = line.split()[4] #extract newick string only (Relate adds some info beforehand)
			ts = from_newick(string,min_edge_length=0.000001) #convert to tskit "tree sequence" (only one tree), it gives errors with branches=0 so I add a minimum b length

			tree = ts.first() #the only tree

        		# shared times
			samples = [int(ts.node(node).metadata['name']) for node in ts.samples()] #get index of each sample in list we gave to relate
			sample_order = np.argsort(samples) #get indices to put in ascending order
			ordered_samples = [ts.samples()[i] for i in sample_order] #order samples as in relate
			sts = get_shared_times(tree, ordered_samples) #get shared times between all pairs of samples, with rows and columns ordered as in relate
			stss.append(sts)

       	 #coalescence times
			cts = sorted([tree.time(i) for i in tree.nodes() if not tree.is_sample(i)]) #coalescence times, in ascending order
			ctss.append(cts)




	np.save(output+".sts.npy", np.array(stss))
	np.save(output+".cts.npy", np.array(ctss))


	stss = np.load(output+".sts.npy")
	_,n,_ = stss.shape



	stss_inv = []
	stss_logdet = []
	smplsss = []
	for sts in stss:
		# chop
		sts_chopped, smpls = chop_shared_times(sts, tCutoff=tCutoff) #shared times and samples of each subtree
		sts_inv = []
		sts_logdet = []
		smplss = []
		# process subtrees
		for st,sm in zip(sts_chopped, smpls):
			stc = center_shared_times(st) #mean center
			stc_inv = np.linalg.pinv(stc) #invert
			stc_logdet = np.linalg.slogdet(st)[1] #log determinant
			sts_inv.append(stc_inv)
			sts_logdet.append(stc_logdet) 
			smplss.append(sm) #samples
		stss_inv.append(sts_inv)
		stss_logdet.append(sts_logdet)
		smplsss.append(smplss)
	# save
	np.save(output+".tcut.sts_mc-inv.npy", np.array(stss_inv, dtype=object))
	np.save(output+".tcut.sts_logdet.npy", np.array(stss_logdet, dtype=object))
	np.save(output+".tcut.samples.npy", np.array(smplsss, dtype=object))


# get variable Nee
	epochs = np.genfromtxt(coal, skip_header=1, skip_footer=1) #time at which each epoch starts (and the final one ends)
	Nes = 1/np.genfromtxt(coal, skip_header=2)[2:] #effective population size during each epoch ==> 1/coal rate for haploids?




	ctss = np.load(output+".cts.npy")


	btss = []
	lpcs = []

	for cts in ctss: 
		# get branching times in ascending order
		T = cts[-1] #TMRCA (to be replaced with tCutoff)
		if tCutoff is not None:
			if tCutoff < T:
				T = tCutoff
		bts = T - np.flip(cts) #branching times, in ascending order
		bts = bts[bts>0] #remove branching times at or before T
		bts = np.append(bts,T) #append total time as last item
		btss.append(bts)
		# get probability of coalescence times under panmictic coalescent with variable Ne
		lpc = _log_coal_density(times=cts, Nes=Nes, epochs=epochs, tCutoff=tCutoff)
		lpcs.append(lpc)
	# save
	np.save(output+".tcut.bts.npy", np.array(btss, dtype=object))
	np.save(output+".tcut.lpcs.npy", np.array(lpcs, dtype=object))


##########################################################################################
##########################################################################################

parser = argparse.ArgumentParser()

#parser.add_argument('-o','--output', metavar='', help='stem name',type=str,nargs=1)
parser.add_argument('-coal','--coal', metavar='', help='.coal file', type =str,nargs=1)
parser.add_argument('-tcut','--tcutoff', metavar='', help='time cutoff for analysis (in generations before present) ',type=float,nargs=1)
parser.add_argument('-data_d','--data_dir', metavar='', help='path to folder with newick trees result of step1 ',type=str,nargs=1)
parser.add_argument('-c','--cores', metavar='', help='number of cores for parallelization', type =int,nargs=1)


#parser.add_argument('-sep_chr' ,'--separate_chr', default= False, help='separate output files by chromosome', action='store_true')
arguments = parser.parse_args()


coal= arguments.coal[0]

dir= arguments.data_dir[0]
search_pattern = os.path.join(dir, "*.newick")
files = glob.glob(search_pattern)
tCutoff=arguments.tcutoff[0]
c=arguments.cores[0]


start_time = time.time()

Parallel(n_jobs=c)(delayed(process_newick)(file,tCutoff,coal) for file in files)
print(f"Time with parallelization: {time.time() - start_time} seconds")


#start_time = time.time()
#for file in files:

#	process_newick(file,tCutoff,coal)
#print(f"Time without parallelization: {time.time() - start_time} seconds")
