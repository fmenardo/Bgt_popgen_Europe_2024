import numpy as np
import argparse
from spacetrees import mle_dispersal, _sds_rho_to_sigma
from tqdm import tqdm
import pandas as pd
import os
import glob


parser = argparse.ArgumentParser()

parser.add_argument('-out','--output', metavar='', help='output name',type=str,nargs=1)
#parser.add_argument('-coal','--coal', metavar='', help='.coal file', type =str,nargs=1)
#parser.add_argument('-tcut','--tcutoff', metavar='', help='time cutoff for analysis (in generations before present) ',type=float,nargs=1)
parser.add_argument('-samples','--samples_file', metavar='', help='path to Relate samples file ',type=str,nargs=1)
parser.add_argument('-meta','--meta', metavar='', help='path to csv with metainfo ',type=str,nargs=1)
parser.add_argument('-ns','--ns', metavar='', help='number of resampled trees to use',type=int,nargs=1)
parser.add_argument('-data_d','--data_dir', metavar='', help='path to folder with newick trees result of step1 ',type=str,nargs=1)



#parser.add_argument('-sep_chr' ,'--separate_chr', default= False, help='separate output files by chromosome', action='store_true')
arguments = parser.parse_args()


dir= arguments.data_dir[0]
search_pattern = os.path.join(dir, "*.newick")
files = glob.glob(search_pattern)
out = arguments.output[0]

samples_file =arguments.samples_file[0]
metadata_file= arguments.meta[0]
df = pd.read_csv(metadata_file)

locations = df[['Sample.Name','Latitude', 'Longitude']]


# Read the sample file
samp = pd.read_csv(samples_file, delimiter=" ", skiprows=1)

# Filter rows based on matching values in "Samples.Names" column and first column
first_column_values = samp.iloc[:, 0]
matching_rows = locations[locations["Sample.Name"].isin(first_column_values)]

matching_rows = matching_rows.set_index("Sample.Name").loc[first_column_values].reset_index()


reordered_locations= matching_rows[['Longitude','Latitude']].to_numpy()

np.save("locations.npy",reordered_locations)

#print(reordered_locations)

lat = locations = df[['Latitude']]
mean_lat =lat.mean()


locations = np.load("locations.npy")


L = len(files)
M = arguments.ns[0]


stss_mc_inv = []
stss_logdet = []
smplss = []
btss = []
lpcss = []

for file in files:
	stem=file.replace(".newick","")
	sts_mc_inv=np.load(stem+".tcut.sts_mc-inv.npy",allow_pickle=True)[:M]
	stss_mc_inv.append(sts_mc_inv)

for file in files:
	stem=file.replace(".newick","")
	sts_logdet=np.load(stem+".tcut.sts_logdet.npy",allow_pickle=True)[:M]
	stss_logdet.append(sts_logdet) 

for file in files:
	stem=file.replace(".newick","")
	smpls=np.load(stem+".tcut.samples.npy",allow_pickle=True)[:M]
	smplss.append(smpls)

for file in files:
	stem=file.replace(".newick","")
	bts=np.load(stem+".tcut.bts.npy",allow_pickle=True)[:M]
	btss.append(bts)

for file in files:
	stem=file.replace(".newick","")
	lpcs=np.load(stem+".tcut.lpcs.npy",allow_pickle=True)[:M]
	lpcss.append(lpcs)




def callbackF(x):
    #  print('{0: 3.6f}   {1: 3.6f}   {2: 3.6f}'.format(x[0], x[1], x[2])) #if important=False
      print('{0: 3.6f}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}'.format(x[0], x[1], x[2], x[3]))


mle = mle_dispersal(locations=locations, shared_times_inverted=stss_mc_inv, log_det_shared_times=stss_logdet, samples=smplss, 
                        sigma0=None, phi0=None, #make educated guess based on first tree at each locus
                        callbackF=callbackF, 
                        important=True, branching_times=btss, logpcoals=lpcss)

print("\n\n")

print('\n',mle)
print('\n',mle[0])
print('\n',mle[1])
print('\n',mle[2])
print('\n',mle[3])




np.save(out+"."+str(M)+"_trees.tcut.mle", mle)


