"""
MODIFIED TO MAP TO SPECIFIC REGION OF GENOME 

The following code takes paired end fastq files and performs quality control (fastp), mapping (bwa), and haplotype-calling (GATK) for each sample. It then sections the GATK produced vcf file by chromosome for further processing. It also subsets the mitochondrion from the final bam file.

Note:
1. Input files are named as *R1(or 2).fastq.gz. 
2. The reference fasta file should have a dictionary and index file (for GATK) and the bwa ref files in the same directory.

"""

import argparse
import os
import subprocess
from datetime import datetime
import re

now = datetime.now()
# dd-mm-YY_H-M-S
today = now.strftime("%d-%m-%Y_%H-%M-%S")

### edit paths as required

#### software ####
fastp_path = "/home/jjigis/bion_tools/fastp"
bwa_path = "/home/jjigis/data/bwa_dir/bwa"
samtools_path = "/home/jjigis/data/samtools-1.17/samtools"
gatk_path = "/home/jjigis/data/gatk-4.4.0.0/gatk"
bcftools_path = "/home/jjigis/data/conda/envs/vcf/bin/bcftools"
tabix_path = "/home/jjigis/data/conda/envs/vcf/bin/tabix"

#### data directories ####
scratch_path = "/home/jjigis/scratch/"
##mito_fol = "/home/jjigis/projects/project_data_prep/analysis/mito_alignments/2023/"
aln_path = "/home/jjigis/projects/project_avrpm17/map_call_pl/bam_files/"
vcf_fol = "/home/jjigis/projects/project_avrpm17/map_call_pl/hap_call_vcf/"
summary_stats_path = "/home/jjigis/projects/project_avrpm17/map_call_pl/summ_stats/"


##os.chdir("../data")

def qc(i1,i2,scratch_fol,minlen,fw,fq,rw,rq): #filenames in the format xyz_1.fastq.gz
    '''takes raw sequences and performs quality trimming & merges overlapping reads'''
    cmd_qc1 = f"{fastp_path} -i {raw_seq_path}{str(i1)} -I {raw_seq_path}{str(i2)} --length_required {str(minlen)} --cut_front --cut_front_window_size {str(fw)} --cut_front_mean_quality {str(fq)} --cut_right --cut_right_window_size {str(rw)} --cut_right_mean_quality {str(rq)} -o {scratch_fol}{str(i1[:-11])}qc1.fq.gz -O {scratch_fol}{str(i1[:-11])}qc2.fq.gz --html {scratch_fol}{str(i1[:-11])}qc_plot.html --json {scratch_fol}{str(i1[:-11])}qc_plot.json"
    subprocess.run(cmd_qc1, shell=True)  
    cmd_qc2 = f"{fastp_path} -i {scratch_fol}{str(i1[:-11])}qc1.fq.gz -I {scratch_fol}{str(i1[:-11])}qc2.fq.gz --merge --overlap_len_require 15 --overlap_diff_limit 1000000 --overlap_diff_percent_limit 10 --merged_out {scratch_fol}{str(i1[:-11])}merged.fq.gz --out1 {scratch_fol}{str(i1[:-11])}out1.fq.gz --out2 {scratch_fol}{str(i2[:-11])}out2.fq.gz --unpaired1 {scratch_fol}{str(i1[:-11])}unpaired1.fq.gz --unpaired2 {scratch_fol}{str(i2[:-11])}unpaired2.fq.gz --json /dev/null --html /dev/null"
    subprocess.run(cmd_qc2, shell=True)  
    merged_wc = subprocess.run(f"zcat {scratch_fol}{str(i1[:-11])}merged.fq.gz | wc -l", shell=True, stdout=subprocess.PIPE)
    return int(merged_wc.stdout)/4   
    
def mapping(ref,i1,scratch_fol):
    '''maps the filtered/trimmed fastq files to reference using bwa'''
    ### number of threads hardcoded. better to pass it as argument 
    cmd_map1 = f"{bwa_path} mem -t 1 {str(ref)} {scratch_fol}{str(i1[:-11])}merged.fq.gz > {scratch_fol}aln_{str(i1[:-11])}merged.sam"
    cmd_map2 = f"{bwa_path} mem -t 1 {str(ref)} {scratch_fol}{str(i1[:-11])}out1.fq.gz {scratch_fol}{str(i1[:-11])}out2.fq.gz > {scratch_fol}aln_{str(i1[:-11])}out1_out2.sam"
    cmd_map3 = f"{bwa_path} mem -t 1 {str(ref)} {scratch_fol}{str(i1[:-11])}unpaired1.fq.gz > {scratch_fol}aln_{str(i1[:-11])}unpaired1.sam"
    cmd_map4 = f"{bwa_path} mem -t 1 {str(ref)} {scratch_fol}{str(i1[:-11])}unpaired2.fq.gz > {scratch_fol}aln_{str(i1[:-11])}unpaired2.sam"
    subprocess.run(cmd_map1, shell=True)
    subprocess.run(cmd_map2, shell=True)
    subprocess.run(cmd_map3, shell=True)
    subprocess.run(cmd_map4, shell=True)

def prop_mapped(bam):
    '''calculates proportion of reads mapped'''
    cmd_pm1 = f"{samtools_path} view -c {str(bam)}" #total number of reads
    cmd_pm2 = f"{samtools_path} view -c -F 260 {str(bam)}" #number of mapped primary aligned reads
    total = subprocess.run(cmd_pm1, shell=True, stdout=subprocess.PIPE)
    mapped = subprocess.run(cmd_pm2, shell=True, stdout=subprocess.PIPE)
    return int(total.stdout), int(mapped.stdout)/int(total.stdout)
    
def index(bam):
    '''creates index of bam file'''
    cmd_i1=f"{samtools_path} index {str(bam)}"
    subprocess.run(cmd_i1, shell=True)
    
def process_bam(i1,scratch_fol):
    '''sorts and formats bam file'''
    cmd_pb1 = f"ls {scratch_fol} | grep sam > {scratch_fol}temp1.txt"
    subprocess.run(cmd_pb1, shell=True)
    with open(f'{scratch_fol}temp1.txt','r') as sam_files:
        lynes = sam_files.readlines()
    sam_files.close()
    for j in lynes:
        sam = f"{scratch_fol}{j.strip()}"
        cmd_pb2 = f"{samtools_path} sort {str(sam)} -o {str(sam[:-3])}bam"
        subprocess.run(cmd_pb2, shell=True)
    cmd_pb3 = f"{samtools_path} merge {scratch_fol}aln_{str(i1[:-11])}compiled.bam {scratch_fol}aln_{str(i1[:-11])}merged.bam {scratch_fol}aln_{str(i1[:-11])}out1_out2.bam {scratch_fol}aln_{str(i1[:-11])}unpaired1.bam {scratch_fol}aln_{str(i1[:-11])}unpaired2.bam"
    subprocess.run(cmd_pb3, shell=True)
    index(f"{scratch_fol}aln_{str(i1[:-11])}compiled.bam")
    total_reads,prop_reads_mapped = prop_mapped(f"{scratch_fol}aln_{str(i1[:-11])}compiled.bam")
    subprocess.run(f"rm {scratch_fol}temp1.txt", shell=True) 
    return total_reads,prop_reads_mapped

def mark_dup(i1, scratch_fol):
    '''first adds fields like read group as required by GATK and then marks duplicates'''
    if "renamed" in str(i1):
        rgsm = str(i1)[8:-12]
    else:
        rgsm = str(i1)[:-12]
    cmd_md1 = f"{gatk_path} --java-options '-Xmx4g' AddOrReplaceReadGroups -I {scratch_fol}aln_{str(i1[:-11])}compiled.bam -O {scratch_fol}aln_{str(i1[:-11])}compiled_RG.bam -RGLB lib_{str(i1[:-12])} -RGPL illumina -RGPU PU_{str(i1[:-12])} -RGSM {rgsm}"
    subprocess.run(cmd_md1, shell = True)
    index(f"{scratch_fol}aln_{str(i1[:-11])}compiled_RG.bam")
    cmd_md2 = f"{gatk_path} --java-options '-Xmx4g' MarkDuplicatesSpark -I {scratch_fol}aln_{str(i1[:-11])}compiled_RG.bam -O {scratch_fol}aln_{str(i1[:-11])}compiled_marked_dup.bam -M {scratch_fol}aln_{str(i1[:-11])}marked_dup_metrics.txt --conf 'spark.local.dir={scratch_fol}'"
    subprocess.run(cmd_md2, shell = True)
    index(f"{scratch_fol}aln_{str(i1[:-11])}compiled_marked_dup.bam")
    
def avg_coverage(scratch_fol, i1):
    '''calculates average genome wide coverage''' 
    cmd_ac1 = f"{samtools_path} coverage -A {scratch_fol}aln_{str(i1[:-11])}compiled_marked_dup.bam -o {summary_stats_path}{str(i1[:-11])}avg_coverage"
    subprocess.run(cmd_ac1, shell=True)
    '''
    cmd_ac2 = f"grep 'Mean coverage' {scratch_fol}{str(i1[:-11])}avg_coverage > {scratch_fol}{str(i1[:-11])}temp_coverage.txt"
    subprocess.run(cmd_ac2, shell=True)
    with open (f"{scratch_fol}{str(i1[:-11])}temp_coverage.txt", "r") as covg_stats:
        covg = covg_stats.readlines()
    covg_stats.close()
    sum_covg = 0
    for i in range(11):
       x_stripped = covg[i].strip()[:-1]
       numm = re.search(r'(\d*\.?\d*)$', x_stripped)
        a = float(numm.group())
        sum_covg += a
    mean_covg = sum_covg/11
    mean_covg = f"{str(round(mean_covg, 2))}x" 
    subprocess.run(f"rm {scratch_fol}{str(i1[:-11])}temp_coverage.txt", shell=True)
    return mean_covg
    '''
    

def hap_caller(ref,scratch_fol,i1):
    '''GATK HaplotypeCaller with default settings'''
    cmd_hc1 = f"{gatk_path} --java-options '-Xmx4g' HaplotypeCaller -R {ref} -I {scratch_fol}aln_{str(i1[:-11])}compiled_marked_dup.bam -O {vcf_fol}{str(i1[:-12])}.avrpm17.bp.vcf.gz -ploidy 1 -ERC BP_RESOLUTION"
    subprocess.run(cmd_hc1, shell = True)
'''        
def section_vcf(i1,scratch_fol):
    chrom_list = ["LR026984.1_chr1","LR026985.1_chr2","LR026986.1_chr3","LR026987.1_chr4","LR026988.1_chr5","LR026989.1_chr6","LR026990.1_chr7","LR026991.1_chr8","LR026992.1_chr9","LR026993.1_chr10","LR026994.1_chr11","MT880591.1","LR026995.1_Un","Bgt_MAT_1_1_3"]
    for chrom in chrom_list:
        cmd_sv1 = f"{bcftools_path} view -r {chrom} {vcf_fol}{str(i1[:-12])}.bp.vcf.gz -Oz -o {vcf_fol}{str(i1[:-11])}{chrom}.bp.vcf.gz"
        subprocess.run(cmd_sv1, shell = True)
        subprocess.run(f"{tabix_path} -p vcf {vcf_fol}{str(i1[:-11])}{chrom}.bp.vcf.gz", shell = True)
    
def extract_mito(i1, scratch_fol):
   # extract mitochondrion from whole genome alignment
    cmd_em1 = f"{samtools_path} view {scratch_fol}aln_{str(i1[:-11])}compiled_marked_dup.bam MT880591.1 -b > {mito_fol}aln_{str(i1[:-11])}mito.bam"
    subprocess.run(cmd_em1, shell = True)
    index(f"{mito_fol}aln_{str(i1[:-11])}mito.bam")
'''    

def extract_mapped(i1,scratch_fol):
    '''extract only mapped reads from bam file'''
    cmd_em1 = f"{samtools_path} view -b -F 4 {scratch_fol}aln_{str(i1[:-11])}compiled_marked_dup.bam > {aln_path}mapped_aln_{str(i1[:-11])}compiled_marked_dup.bam"
    subprocess.run(cmd_em1, shell = True)
    index(f"{aln_path}mapped_aln_{str(i1[:-11])}compiled_marked_dup.bam")

#### argument parser    
parser = argparse.ArgumentParser()

parser.add_argument('-ref', '--reference', metavar='', help='absolute path to the reference fasta file')
parser.add_argument('-minlen','--minimum_read_length', type=int, default='', metavar = '', help='length_required in fastp', nargs=1)
parser.add_argument('-rw','--right_window_fastp', type=int, default='', metavar = '', help='cut_right_window_size in fastp', nargs=1)
parser.add_argument('-fw','--front_window_fastp', type=int, default='', metavar = '', help='cut_front_window_size in fastp', nargs=1)
parser.add_argument('-rq','--right_mean_qual_fastp', type=int, default='', metavar='',  help='cut_right_mean_quality in fastp', nargs=1)
parser.add_argument('-fq','--front_mean_qual_fastp', type=int, default='', metavar='',  help='cut_front_mean_quality in fastp', nargs=1)
parser.add_argument('-i', '--input', metavar = '', help='absolute path to input fastq file. Do not use "~" notation')

ag = parser.parse_args()

infile_path = re.split('\/',ag.input)
fastq_file = infile_path[-1]
i1 = f"{str(fastq_file)[:-11]}R1.fastq.gz"
i2 = f"{str(fastq_file)[:-11]}R2.fastq.gz"  
raw_seq_path = f"{str('/'.join(infile_path[:-1]))}/"

ref = str(ag.reference)

fl = f"{i1[:-12]}/{today}/"
scratch_folder= os.path.join(scratch_path,fl)
os.makedirs(scratch_folder)

    
merged = qc(i1,i2,scratch_folder,ag.minimum_read_length[0],ag.front_window_fastp[0],ag.front_mean_qual_fastp[0],ag.right_window_fastp[0],ag.right_mean_qual_fastp[0])
mapping(ref,i1,scratch_folder)
total_reads, prop_mapped = process_bam(i1,scratch_folder)
mark_dup(i1,scratch_folder)
avg_coverage(scratch_folder,i1)
hap_caller(ref,scratch_folder,i1)
extract_mapped(i1,scratch_folder)
#section_vcf(i1,scratch_folder)
#extract_mito(i1, scratch_folder)

totreads = re.findall(r'\d+',str(total_reads))[0]

with open(f"{summary_stats_path}{str(i1[:-12])}_pl_stats.csv", 'w') as summ_stats:
    summ_stats.write(f"isolate,{i1[:-12]}\ntotal number of reads,{totreads}\nnumber of reads merged due to overlap,{merged}\nproportion of reads mapped,{str(prop_mapped)}\n")
#average genome wide coverage,{str(mean_coverage).strip()[:-1]}\n")

now = datetime.now()
end_time = now.strftime("%d-%m-%Y_%H-%M-%S")

print(f"start time = {today}")
print(f"end time = {end_time}")

