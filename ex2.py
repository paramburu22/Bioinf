# Ex2.py
from Bio.Blast import NCBIWWW
from Bio import SeqIO

def main():
    input_file = "HBB_correct_ORF.fasta"
    output_file = "blast_HBB_remote.xml"

    print(f"🔍 Leyendo secuencia desde {input_file} ...")
    record = SeqIO.read(input_file, format="fasta")

    print("🚀 Ejecutando BLASTP remoto contra base de datos 'swissprot' ...")
    result_handle = NCBIWWW.qblast(
        program="blastp",       # tipo de BLAST
        database="swissprot",   # base de datos de proteínas
        sequence=record.seq,    # secuencia de entrada
        hitlist_size=10         # limitar resultados
    )

    print(f"💾 Guardando resultados en {output_file} ...")
    with open(output_file, "w") as out_f:
        out_f.write(result_handle.read())
    result_handle.close()

    print("✅ BLAST remoto completado correctamente.")

if __name__ == "__main__":
    main()
