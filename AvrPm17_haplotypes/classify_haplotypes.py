from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def classify_haps(in_file, start, end):
    sequence = []
    haps = {}
    for record in SeqIO.parse(in_file, "fasta"):
        my_seq = record.seq[start -1 : end]
        if my_seq in sequence:
            hapid = sequence.index(my_seq)
            haps[hapid].append(record.id)
        else:
            sequence.append(my_seq)
            hapid = sequence.index(my_seq)
            haps[hapid] = [record.id]
    return haps



## classify RC-DNA of whole coding sequence
whole_cds_RC = classify_haps("avrpm17_ext_eur_only_cds_RC_unique_no_miss.fa", 1, 330)
with open ("hap_class_cds_RC_no_miss.csv", "w+") as f:
    f.write("hap_class,hap_id\n")
    for hap in whole_cds_RC.keys():
        for ids in whole_cds_RC[hap]:
            f.write(f"{hap},{ids}\n")
            

## classify protein sequence of whole cds
whole_cds_prot = classify_haps("avrpm17_ext_eur_only_cds_trns_unique_no_miss.fa", 1, 110)
with open ("hap_class_cds_trns_no_miss.csv", "w+") as f:
    f.write("hap_class,hap_id\n")
    for hap in whole_cds_prot.keys():
        for ids in whole_cds_prot[hap]:
            f.write(f"{hap},{ids}\n")


## classify RC-DNA in mature protein seq only
mp_RC = classify_haps("avrpm17_ext_eur_only_cds_RC_unique_no_miss.fa", 73, 330)
with open ("hap_class_mp_RC_no_miss.csv", "w+") as f:
    f.write("hap_class,hap_id\n")
    for hap in mp_RC.keys():
        for ids in mp_RC[hap]:
            f.write(f"{hap},{ids}\n")

## classify protein sequence of mature protein only
mp_prot = classify_haps("avrpm17_ext_eur_only_cds_trns_unique_no_miss.fa", 25, 110)
with open ("hap_class_mp_trns_no_miss.csv", "w+") as f:
    f.write("hap_class,hap_id\n")
    for hap in mp_prot.keys():
        for ids in mp_prot[hap]:
            f.write(f"{hap},{ids}\n")
