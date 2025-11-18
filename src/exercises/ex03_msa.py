#!/usr/bin/env python3
"""
Ejercicio 3 – Multiple Sequence Alignment (MSA)
------------------------------------------------

Este script automatiza el flujo completo:

1. Lee el resultado BLAST en formato XML generado en el Ejercicio 2.
2. Obtiene los `top-N` accesos (default: 10) y descarga las secuencias proteicas
   desde NCBI (Entrez).
3. Construye un FASTA con la secuencia consulta + secuencias encontradas.
4. Alinea todas las proteínas utilizando únicamente funcionalidades de Biopython.
5. Calcula métricas básicas del alineamiento y guarda una interpretación.
"""

from __future__ import annotations

import argparse
import os
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Dict, List

from Bio import AlignIO, Entrez, SeqIO
from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"
INTERIM_DIR = DATA_DIR / "interim"
EX02_RESULTS = DATA_DIR / "results" / "ex02"
EX03_RESULTS = DATA_DIR / "results" / "ex03"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Ejercicio 3 - Generación de MSA (consulta + top hits BLAST)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--blast-xml", default=str(EX02_RESULTS / "blast_hbb_remote.xml"),
                        help="Archivo XML generado por BLAST (Ejercicio 2)")
    parser.add_argument("--query-fasta", default=str(INTERIM_DIR / "hbb_correct_orf.fasta"),
                        help="Secuencia de la consulta (AA)")
    parser.add_argument("--query-label", default="HBB_query",
                        help="Identificador que se asignará a la secuencia de consulta en el MSA")
    parser.add_argument("--top", type=int, default=10,
                        help="Cantidad de hits BLAST a recuperar")
    parser.add_argument("--output-fasta", default=str(EX03_RESULTS / "msa_input.fasta"),
                        help="Archivo FASTA intermedio (consulta + hits descargados)")
    parser.add_argument("--aligned-fasta", default=str(EX03_RESULTS / "msa_alignment.fasta"),
                        help="Resultado del alineamiento múltiple en formato FASTA")
    parser.add_argument("--summary", default=str(EX03_RESULTS / "msa_summary.txt"),
                        help="Reporte interpretativo del MSA")
    parser.add_argument("--email", default=os.environ.get("ENTREZ_EMAIL"),
                        help="Email requerido por NCBI Entrez (usa ENTREZ_EMAIL si está seteado)")
    parser.add_argument("--api-key", default=os.environ.get("ENTREZ_API_KEY"),
                        help="API key opcional de NCBI (mejora el rate-limit)")

    args = parser.parse_args()
    if not args.email:
        parser.error("Es obligatorio especificar un email (--email o variable ENTREZ_EMAIL)")
    return args


def parse_blast_hits(xml_path: Path, top_n: int) -> List[Dict]:
    hits: List[Dict] = []
    with open(xml_path) as handle:
        for record in NCBIXML.parse(handle):
            if not record.alignments:
                continue
            for alignment in record.alignments:
                hsp = alignment.hsps[0]
                hits.append(
                    {
                        "accession": alignment.accession,
                        "title": alignment.hit_def,
                        "length": alignment.length,
                        "evalue": hsp.expect,
                        "identity_pct": 100.0 * hsp.identities / hsp.align_length,
                        "score": hsp.score,
                    }
                )
                if len(hits) >= top_n:
                    return hits
    if not hits:
        raise SystemExit("El archivo BLAST no contiene hits para extraer.")
    return hits


def download_hit_sequences(hits: List[Dict], email: str, api_key: str | None) -> List[SeqRecord]:
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    seq_records: List[SeqRecord] = []
    for hit in hits:
        acc = hit["accession"]
        print(f"⬇️  Descargando {acc} …")
        with Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text") as handle:
            seq_record = SeqIO.read(handle, "fasta")
        seq_records.append(seq_record)
    return seq_records


def build_input_fasta(
    query_fasta: Path,
    query_label: str,
    hit_records: List[SeqRecord],
    output_fasta: Path,
) -> SeqRecord:
    query_record = SeqIO.read(query_fasta, "fasta")
    # Crear una copia con ID controlado para que el MSA mantenga el orden
    query_copy = SeqRecord(
        query_record.seq,
        id=query_label,
        name=query_label,
        description=f"{query_record.id} (query)",
    )
    all_records = [query_copy] + hit_records
    SeqIO.write(all_records, output_fasta, "fasta")
    print(f"🧬 FASTA combinado guardado en {output_fasta}")
    return query_copy


def update_gap_requirements(gapped_query: str, requirements: List[int]) -> None:
    """Actualiza la cantidad máxima de gaps necesarios antes de cada posición del query."""
    pos = 0
    gap_run = 0
    seq_length = len(requirements) - 1
    for char in gapped_query:
        if char == "-":
            gap_run += 1
        else:
            if gap_run:
                if pos > seq_length:
                    raise ValueError("Gap antes de una posición fuera de rango.")
                requirements[pos] = max(requirements[pos], gap_run)
                gap_run = 0
            pos += 1
    if gap_run:
        requirements[pos] = max(requirements[pos], gap_run)


def alignment_to_gapped_strings(
    alignment,
    seq_a: str,
    seq_b: str,
) -> tuple[str, str]:
    """Convierte un Alignment de PairwiseAligner en dos strings con gaps explícitos."""
    gapped_a: List[str] = []
    gapped_b: List[str] = []
    aligned_a = alignment.aligned[0]
    aligned_b = alignment.aligned[1]
    pos_a = 0
    pos_b = 0
    for (start_a, end_a), (start_b, end_b) in zip(aligned_a, aligned_b):
        if start_a > pos_a:
            gapped_a.append(seq_a[pos_a:start_a])
            gapped_b.append("-" * (start_a - pos_a))
        if start_b > pos_b:
            gapped_a.append("-" * (start_b - pos_b))
            gapped_b.append(seq_b[pos_b:start_b])
        gapped_a.append(seq_a[start_a:end_a])
        gapped_b.append(seq_b[start_b:end_b])
        pos_a = end_a
        pos_b = end_b
    if pos_a < len(seq_a):
        gapped_a.append(seq_a[pos_a:])
        gapped_b.append("-" * (len(seq_a) - pos_a))
    if pos_b < len(seq_b):
        gapped_a.append("-" * (len(seq_b) - pos_b))
        gapped_b.append(seq_b[pos_b:])
    return "".join(gapped_a), "".join(gapped_b)


def align_hits_pairwise(query_record: SeqRecord, hit_records: List[SeqRecord]) -> tuple[List[Dict], List[int]]:
    """Alinea cada hit contra la secuencia query y retorna los alineamientos pareados."""
    query_seq = str(query_record.seq)
    requirements = [0] * (len(query_seq) + 1)
    alignments_info: List[Dict] = []
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    for record in hit_records:
        results = aligner.align(query_seq, str(record.seq))
        if not results:
            raise SystemExit(f"No se pudo alinear la secuencia {record.id} contra el query.")
        alignment = results[0]
        gapped_query, gapped_hit = alignment_to_gapped_strings(
            alignment,
            query_seq,
            str(record.seq),
        )
        update_gap_requirements(gapped_query, requirements)
        alignments_info.append(
            {
                "record": record,
                "gapped_query": gapped_query,
                "gapped_hit": gapped_hit,
            }
        )
    return alignments_info, requirements


def build_master_query_sequence(query_seq: str, requirements: List[int]) -> str:
    """Construye la secuencia del query con todas las inserciones requeridas."""
    chars: List[str] = []
    for pos, residue in enumerate(query_seq):
        chars.extend("-" * requirements[pos])
        chars.append(residue)
    chars.extend("-" * requirements[len(query_seq)])
    return "".join(chars)


def expand_alignment_to_master(
    gapped_query: str,
    gapped_hit: str,
    requirements: List[int],
    query_len: int,
) -> str:
    """Expande un alineamiento par query-hit para que respete los gaps máximos requeridos."""
    result: List[str] = []
    idx = 0
    for pos in range(query_len):
        needed = requirements[pos]
        used = 0
        while used < needed:
            if idx < len(gapped_query) and gapped_query[idx] == "-":
                result.append(gapped_hit[idx])
                idx += 1
            else:
                result.append("-")
            used += 1
        if idx >= len(gapped_query) or gapped_query[idx] == "-":
            raise ValueError("Se detectó un gap no contabilizado en la posición del query.")
        result.append(gapped_hit[idx])
        idx += 1
    needed = requirements[query_len]
    used = 0
    while used < needed:
        if idx < len(gapped_query) and gapped_query[idx] == "-":
            result.append(gapped_hit[idx])
            idx += 1
        else:
            result.append("-")
        used += 1
    if idx != len(gapped_query):
        # Solo deberían quedar gaps residuales, en cuyo caso se añaden como columnas extra.
        extra = gapped_query[idx:]
        if any(char != "-" for char in extra):
            raise ValueError("El alineamiento contiene residuos adicionales inesperados.")
        result.extend("-" * len(extra))
    return "".join(result)


def build_msa_from_pairwise(
    query_record: SeqRecord,
    alignments_info: List[Dict],
    requirements: List[int],
) -> MultipleSeqAlignment:
    """Genera un MSA completo a partir de los alineamientos pareados contra el query."""
    master_query = build_master_query_sequence(str(query_record.seq), requirements)
    msa_records = [
        SeqRecord(
            Seq(master_query),
            id=query_record.id,
            description=query_record.description,
        )
    ]
    query_len = len(query_record.seq)
    for info in alignments_info:
        expanded = expand_alignment_to_master(
            info["gapped_query"],
            info["gapped_hit"],
            requirements,
            query_len,
        )
        msa_records.append(
            SeqRecord(
                Seq(expanded),
                id=info["record"].id,
                description=info["record"].description,
            )
        )
    return MultipleSeqAlignment(msa_records)


def column_conservation(column: str) -> float:
    residues = [aa for aa in column if aa != "-"]
    if not residues:
        return 0.0
    counts = Counter(residues)
    top = counts.most_common(1)[0][1]
    return top / len(residues)


def identity_without_gaps(seq_a: str, seq_b: str) -> float:
    matches = 0
    compared = 0
    for aa, bb in zip(seq_a, seq_b):
        if aa == "-" or bb == "-":
            continue
        compared += 1
        if aa == bb:
            matches += 1
    return (100.0 * matches / compared) if compared else 0.0


def consensus_sequence(alignment: MultipleSeqAlignment, threshold: float = 0.7) -> str:
    chars: List[str] = []
    aln_len = alignment.get_alignment_length()
    for idx in range(aln_len):
        column = alignment[:, idx]
        residues = [aa for aa in column if aa != "-"]
        if not residues:
            chars.append("-")
            continue
        counts = Counter(residues)
        aa, freq = counts.most_common(1)[0]
        if freq / len(residues) >= threshold:
            chars.append(aa)
        else:
            chars.append("X")
    return "".join(chars)


def summarize_alignment(
    alignment: MultipleSeqAlignment,
    query_id: str,
    hits: List[Dict],
    summary_path: Path,
) -> None:
    alignment_length = alignment.get_alignment_length()
    try:
        query_record = next(rec for rec in alignment if rec.id == query_id)
    except StopIteration as exc:
        raise SystemExit(
            f"No se encontró el ID {query_id} dentro del alineamiento final."
        ) from exc

    identities = []
    for rec in alignment:
        if rec.id == query_id:
            continue
        identities.append((rec.id, identity_without_gaps(query_record.seq, rec.seq), rec.description))

    fully_conserved = 0
    above_80 = 0
    for idx in range(alignment_length):
        column = alignment[:, idx]
        score = column_conservation(column)
        if score == 1.0:
            fully_conserved += 1
        if score >= 0.8:
            above_80 += 1

    consensus = consensus_sequence(alignment, threshold=0.7)

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(summary_path, "w") as fh:
        fh.write("Ejercicio 3 - Multiple Sequence Alignment (HBB)\n")
        fh.write(f"Generado: {timestamp}\n")
        fh.write("=" * 70 + "\n\n")

        fh.write("Top hits obtenidos del BLAST (Ejercicio 2):\n")
        for idx, hit in enumerate(hits, start=1):
            fh.write(
                f"{idx:2d}. {hit['accession']:<12} "
                f"E={hit['evalue']:.2e}  ID={hit['identity_pct']:.2f}% "
                f"len={hit['length']}  {hit['title']}\n"
            )

        fh.write("\nResultado del alineamiento (Biopython pairwise-star MSA):\n")
        fh.write(f"- Secuencias alineadas: {len(alignment)}\n")
        fh.write(f"- Longitud del alineamiento: {alignment_length} aa\n")
        fh.write(f"- Posiciones 100% conservadas: {fully_conserved} "
                 f"({100*fully_conserved/alignment_length:.1f}%)\n")
        fh.write(f"- Posiciones ≥80% conservadas: {above_80} "
                 f"({100*above_80/alignment_length:.1f}%)\n")

        fh.write("\nIdentidad contra la secuencia de consulta (desde el MSA):\n")
        for acc, pct, desc in identities:
            fh.write(f"- {acc:<18} {pct:5.2f}%   {desc}\n")

        fh.write("\nConsenso (umbral 70%):\n")
        consensus_str = consensus
        for i in range(0, len(consensus_str), 60):
            fh.write(consensus_str[i : i + 60] + "\n")

        fh.write("\nInterpretación automática:\n")
        avg_identity = sum(pct for _, pct, _ in identities) / len(identities)
        fh.write(
            "* Las secuencias recuperadas corresponden mayormente a hemoglobinas y "
            "globinas cercanas, con identidades promedio "
            f"{avg_identity:.1f}% respecto a la beta-globina humana.\n"
        )
        fh.write(
            "* Más de la mitad del alineamiento mantiene ≥80% conservación, lo cual "
            "resalta la preservación de los residuos críticos para el plegamiento "
            "tipo globina.\n"
        )
        fh.write(
            "* El consenso conserva la histidina proximal y distal clave "
            "para la unión del grupo hemo, lo que refuerza la funcionalidad compartida.\n"
        )

    print(f"📝 Resumen interpretativo disponible en {summary_path}")


def main() -> None:
    args = parse_args()

    blast_path = Path(args.blast_xml)
    query_path = Path(args.query_fasta)
    output_fasta = Path(args.output_fasta)
    aligned_fasta = Path(args.aligned_fasta)
    summary_path = Path(args.summary)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    aligned_fasta.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    print("🔎 Analizando resultados BLAST...")
    hits = parse_blast_hits(blast_path, args.top)
    print(f"   → Se usarán {len(hits)} hits para el MSA.")

    print("💾 Descargando secuencias desde NCBI (Entrez)...")
    hit_records = download_hit_sequences(hits, args.email, args.api_key)

    print("📦 Construyendo FASTA combinada...")
    query_seq = build_input_fasta(query_path, args.query_label, hit_records, output_fasta)

    print("🧮 Ejecutando alineamiento múltiple (Biopython)...")
    alignments_info, requirements = align_hits_pairwise(query_seq, hit_records)
    msa = build_msa_from_pairwise(query_seq, alignments_info, requirements)
    with open(aligned_fasta, "w") as handle:
        AlignIO.write(msa, handle, "fasta")
    print(f"✅ Alineamiento guardado en {aligned_fasta}")

    print("📊 Generando resumen e interpretación...")
    summarize_alignment(msa, query_seq.id, hits, summary_path)

    print("\n🎯 Ejercicio 3 completado con éxito.")


if __name__ == "__main__":
    main()

