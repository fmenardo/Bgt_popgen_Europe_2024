from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def extract_rc_seq(fasta_file, output_file, rc_file, start, end):
    with open(output_file, "a+") as f, open(rc_file, "a+") as rc:
        for record in SeqIO.parse(str(fasta_file), "fasta"):
            sequence = record.seq[start -1 : end]
            seq_record = SeqRecord(sequence, id=record.id, description="")
            SeqIO.write(seq_record, f, "fasta")
            seq_rc_record = SeqRecord(sequence.reverse_complement(), id=f"{record.id}_RC", description = "")
            SeqIO.write(seq_rc_record, rc, "fasta")
            


fasta_file = "avrpm17_msa_ext_eur_no_missing.fa"
output_file = "avrpm17_msa_ext_eur_only_gene_no_missing.fa"
rc_file = "avrpm17_msa_ext_eur_only_gene_RC_no_missing.fa"
start = 2001
end = 2386
            
extract_rc_seq(fasta_file, output_file, rc_file, start, end)
