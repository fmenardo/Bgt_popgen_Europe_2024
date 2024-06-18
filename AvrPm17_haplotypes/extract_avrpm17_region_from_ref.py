from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def extract_sequence(fasta_file, contig_name, start, end):

    sequence = ""
    for record in SeqIO.parse(str(fasta_file), "fasta"):
        if record.id == contig_name:
            sequence = record.seq[start - 1:end]
            break
    return str(sequence)


def write_to_fasta(output_file, sequence, header):

    with open(output_file, "w") as f:
        seq_record = SeqRecord(sequence, id=header, description="")
        SeqIO.write(seq_record, f, "fasta")


fasta_file = "../vcf_project_tritici/GCA_900519115.1_2022_bgt_ref_mating_type.fa"
contig_name = "LR026984.1_chr1"
start= 4365017-2000
end = 4365402+2000

for record in SeqIO.parse(str(fasta_file), "fasta"):
        if record.id == contig_name:
                sequence = record.seq[start - 1:end]
                break

#sequence = extract_sequence(fasta_file, contig_name, start_position, end_position)


output_file = "avrpm17_locus+-2Kb.fa"



write_to_fasta(output_file, sequence,"avrpm17_locus+-2Kb")
