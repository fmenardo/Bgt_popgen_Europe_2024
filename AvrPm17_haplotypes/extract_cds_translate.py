from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

with open("avrpm17_ext_eur_only_cds_no_miss.fa", "w+") as cds, open ("avrpm17_ext_eur_only_cds_RC_no_miss.fa", "w+") as cds_rc, open ("avrpm17_ext_eur_only_cds_trns_no_miss.fa", "w+") as cds_tr:
    for record in SeqIO.parse("avrpm17_msa_ext_eur_only_gene_no_missing.fa", "fasta"):
        sequence = record.seq[1-1 : 95]+record.seq[152-1 : 386]
        seq_record = SeqRecord(sequence, id=f"{record.id}_cds", description="")
        SeqIO.write(seq_record, cds, "fasta")
        rc = sequence.reverse_complement()
        seq_rc_record = SeqRecord(rc, id=f"{record.id}_cds_RC", description = "")
        SeqIO.write(seq_rc_record, cds_rc, "fasta")
        pr = rc.translate()
        seq_pr_record = SeqRecord(pr, id = f"{record.id}_cds_aa", description = "")
        SeqIO.write(seq_pr_record, cds_tr, "fasta")
