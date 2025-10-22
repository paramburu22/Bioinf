# Ex2_remote.py
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

INPUT = "HBB_ORFs.fasta"
OUT_XML = "blast_frames_remote.xml"
SUMMARY = "blast_frames_summary.tsv"

def run_blast_all_frames():
    seqs = list(SeqIO.parse(INPUT, "fasta"))
    # concatenar consultas en un solo XML (una por registro)
    with open(OUT_XML, "w") as xml_out, open(SUMMARY, "w") as summ:
        summ.write("frame\tbest_hit\tevalue\tpct_identity\talign_len\tquery_len\n")
        for rec in seqs:
            print(f"🚀 BLASTP remoto → {rec.id}")
            handle = NCBIWWW.qblast(
                program="blastp",
                database="swissprot",
                sequence=str(rec.seq),
                hitlist_size=10
            )
            xml_text = handle.read()
            xml_out.write(xml_text)
            handle.close()

            # parsear inmediatamente el top hit de este frame
            blast_record = next(NCBIXML.read(open(OUT_XML)))  # lee el último desde archivo
            if blast_record.alignments:
                aln = blast_record.alignments[0]
                hsp = aln.hsps[0]
                pct_id = 100.0 * hsp.identities / hsp.align_length
                summ.write(f"{rec.id}\t{aln.hit_def}\t{hsp.expect:.2e}\t{pct_id:.2f}\t{hsp.align_length}\t{blast_record.query_length}\n")
            else:
                summ.write(f"{rec.id}\t(no hits)\tNA\tNA\t0\t{blast_record.query_length}\n")

if __name__ == "__main__":
    run_blast_all_frames()
    print(f"✅ Listo: {OUT_XML} y {SUMMARY}")
    