import argparse
import numpy as np




parser = argparse.ArgumentParser()

parser.add_argument('-vcf','--vcf_file', metavar='', help='input vcf with no missing data (for fs)', type =str,nargs=1)
parser.add_argument('-o','--output', metavar='', help='output stem',type=str,nargs=1)
#parser.add_argument('-sep_chr' ,'--separate_chr', default= False, help='separate output files by chromosome', action='store_true')
arguments = parser.parse_args()






count=1
# read vcf
matrix=[]
list_positions=[]
with open(arguments.vcf_file[0], 'r') as file:
	for line in file:
		if line.startswith("##"):
			continue
		line_stripped = line.rstrip('\n')
		columns = line_stripped.split('\t')
		row=[]
		indices_to_remove = [2,3,4,5,6,7,8]

		indices_to_remove.sort(reverse=True)

		for index in indices_to_remove:
			columns.pop(index)
		if count==1:
			matrix.append(columns)
			id_list=columns[1]
		else:
			index_of_chr = columns[0].find("chr")
			chr = columns[0][index_of_chr + len("chr"):]
			list_positions.append(chr+","+columns[1])
			for i in range(2, len(columns)):
				columns[i]= columns[i][:1]
			matrix.append(columns)
		count=count+1




matrix_np = np.array(matrix)
list_id=matrix_np[0,2:]
matrix_np_t=np.transpose(matrix_np)

matrix_res=matrix_np_t[2:, 1:]
matrix_np_t[1,0] = "P"


num_rows, num_columns = matrix_res.shape

n_hap= num_rows
n_SNPs= num_columns
pos_line= matrix_np_t[1, :]



file_path = arguments.output[0]+".hap_file"

with open(file_path, 'w') as file:
    file.write(str(n_hap) + '\n')

    file.write(str(n_SNPs) + '\n')

    file.write(' '.join(map(str, pos_line)) + '\n')

    for row in matrix_res:
        file.write(''.join(map(str, row)) + '\n')


file_path = arguments.output[0]+".id_file"

with open(file_path, 'w') as file:

	file.write('\n'.join(map(str, list_id)))


file_path = arguments.output[0]+".pos_file"

with open(file_path, 'w') as file:

        file.write('\n'.join(map(str, list_positions)))
