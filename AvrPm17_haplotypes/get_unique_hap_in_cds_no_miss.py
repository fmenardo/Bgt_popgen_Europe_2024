from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

seq_dict = {}
for record in SeqIO.parse("avrpm17_ext_eur_only_cds_no_miss.fa", "fasta"):
    seq_dict[record.id] = record.seq

with open ("../map_call_pl/het_sample_list_no_miss" , 'r') as sample_file:
    samples = sample_file.readlines()
    for s in samples:
        sample = s.strip()
        if seq_dict[f"{sample}_1_avrpm17_cds"] == seq_dict[f"{sample}_2_avrpm17_cds"]:
            seq_dict.pop(f"{sample}_2_avrpm17_cds")
            
with open ("avrpm17_ext_eur_only_cds_unique_no_miss.fa", "w+") as cds, open("avrpm17_ext_eur_only_cds_RC_unique_no_miss.fa", "w+") as cds_rc, open ("avrpm17_ext_eur_only_cds_trns_unique_no_miss.fa", "w+") as cds_tr:
    for header, seq in seq_dict.items():
        seq_rec = SeqRecord(seq, id=header, description = "")
        SeqIO.write(seq_rec, cds, "fasta")
        rc = seq.reverse_complement()
        seq_rc_record = SeqRecord(rc, id=header, description = "")
        SeqIO.write(seq_rc_record, cds_rc, "fasta")
        pr = rc.translate()
        seq_pr_record = SeqRecord(pr, id=header, description = "")
        SeqIO.write(seq_pr_record, cds_tr, "fasta")
