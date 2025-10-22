# download_blast_hits.py
from Bio.Blast import NCBIXML
from Bio import Entrez, SeqIO

# === CONFIGURACIÓN ===
Entrez.email = "tu_email@ejemplo.com"  # Requerido por NCBI
blast_xml = "blast_HBB_remote.xml"
output_fasta = "MSA_input.fasta"
top_n = 10

# === LEER RESULTADOS DEL BLAST ===
print("🔍 Leyendo resultados de BLAST...")
with open(blast_xml) as result_handle:
    blast_records = NCBIXML.parse(result_handle)
    record = next(blast_records)  # solo el primer query

# === EXTRAER LOS 10 MEJORES HITS ===
hits = []
for alignment in record.alignments[:top_n]:
    acc = alignment.accession
    hits.append(acc)
    print("✅ Hit encontrado:", acc)

# === DESCARGAR SECUENCIAS EN FORMATO FASTA ===
print("\n⬇️ Descargando secuencias desde NCBI...")
with open(output_fasta, "w") as out_f:
    for acc in hits:
        handle = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
        seq_record = handle.read()
        out_f.write(seq_record)
        handle.close()
        print(f"Guardado {acc}")

print(f"\n✅ Archivo final creado: {output_fasta}")
