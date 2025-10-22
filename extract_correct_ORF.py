# Ex1_from_cds.py
from Bio import SeqIO

gbk = "HBB_NM000518.gb"   # tu GenBank
out_fa = "HBB_correct_ORF2.fasta"

record = SeqIO.read(gbk, "genbank")

# Buscar el feature CDS
cds_feats = [f for f in record.features if f.type == "CDS"]
assert cds_feats, "No se encontró feature CDS en el GenBank"

cds = cds_feats[0]
cds_seq = cds.extract(record.seq)                    # nucleótidos de CDS
prot = cds_seq.translate(to_stop=True)               # proteína madura

with open(out_fa, "w") as f:
    f.write(">HBB_protein_correct_ORF\n")
    f.write(str(prot) + "\n")

print("Listo →", out_fa, "len=", len(prot))
