import numpy as np
import argparse
from spacetrees import _sds_rho_to_sigma, _log_birth_density, locate_ancestors
from utils import chop_shared_times
import pandas as pd
import os
import glob
import time
from joblib import Parallel, delayed


def run_locate(file,mle,M,tCutoff,sigma):

	stss_chopped = []
	samples = []

	stem=file.replace(".newick","")
	for sts in np.load(stem+".sts.npy"):
		sts_chopped, smpls = chop_shared_times(sts, tCutoff=tCutoff) #shared times and sample indices in each subtree (here just 1 subtree per tree)
		stss_chopped.append(sts_chopped)
		samples.append(smpls)


	# importance weights
	btss = np.load(stem+".tcut.bts.npy", allow_pickle=True) #birth times
	phi = mle[-1] #mle branching rate
	lbds = np.array([_log_birth_density(bts, phi, n) for bts in btss]) #log probability densities of birth times
	lpcs = np.load(stem+".tcut.lpcs.npy", allow_pickle=True) #log probability densities of coalescence times
	log_weights = lbds - lpcs #log importance weights

	# locate 
	ancestor_locations = locate_ancestors(ancestor_samples=ancestor_samples, ancestor_times=ancestor_times, shared_times_chopped=stss_chopped[:M], samples=samples[:M], locations=locations, log_weights= log_weights[:M],sigma=sigma)

	np.save(stem+".anc-locs-full-set", ancestor_locations)



parser = argparse.ArgumentParser()

#parser.add_argument('-out','--output', metavar='', help='output name',type=str,nargs=1)
#parser.add_argument('-coal','--coal', metavar='', help='.coal file', type =str,nargs=1)
parser.add_argument('-tcut','--tcutoff', metavar='', help='time cutoff for analysis (in generations before present) ',type=float,nargs=1)
#parser.add_argument('-samples','--samples_file', metavar='', help='path to Relate samples file ',type=str,nargs=1)
parser.add_argument('-mle','--mle', metavar='', help='mle file ',type=str,nargs=1)
parser.add_argument('-ns','--ns', metavar='', help='number of times branches are resampled by relate ',type=int,nargs=1)
parser.add_argument('-data_d','--data_dir', metavar='', help='path to folder with newick trees result of step1 ',type=str,nargs=1)
parser.add_argument('-c','--cores', metavar='', help='number of cores for parallelization', type =int,nargs=1)



#parser.add_argument('-sep_chr' ,'--separate_chr', default= False, help='separate output files by chromosome', action='store_true')
arguments = parser.parse_args()


dir= arguments.data_dir[0]
search_pattern = os.path.join(dir, "*.newick")
files = glob.glob(search_pattern)
mle_file = arguments.mle[0]


locations = np.load("locations.npy")
n= len(locations)

ancestor_samples = range(n) #which samples to find ancestors of
tCutoff=float(arguments.tcutoff[0])   #determines which dispersal rate we use and how far back we locate ancestors
max_time = tCutoff
ancestor_times = np.linspace(10,max_time,10) #times to find ancestors
ancestor_times= [1,2,3,4,5,10,25,50,100]

M = arguments.ns[0]
c = arguments.cores[0]

#dispersal rate
mle = np.load(mle_file) #mle dispersal rate and branching rate
sigma = _sds_rho_to_sigma(mle[:-1]) #mle[0], mle[1], mle[2]) #as covariance matrix

start_time = time.time()

Parallel(n_jobs=c)(delayed(run_locate)(file,mle,M,tCutoff,sigma) for file in files)
print(f"Time with parallelization: {time.time() - start_time} seconds")
