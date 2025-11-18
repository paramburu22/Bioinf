#!/usr/bin/env python3
"""Ejercicio 2 - BLAST remoto contra SwissProt."""

import argparse
import io
from pathlib import Path

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_INTERIM = PROJECT_ROOT / "data" / "interim"
RESULTS_DIR = PROJECT_ROOT / "data" / "results" / "ex02"
DEFAULT_INPUT = DATA_INTERIM / "hbb_orfs.fasta"
DEFAULT_XML = RESULTS_DIR / "blast_hbb_remote.xml"
DEFAULT_SUMMARY = RESULTS_DIR / "blast_summary.tsv"
DEFAULT_TEXT = RESULTS_DIR / "blast.out"


def run_blast_all_frames(
    input_fasta: Path,
    xml_path: Path,
    summary_path: Path,
    text_path: Path,
    hitlist_size: int,
) -> None:
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    if not sequences:
        raise SystemExit(f"No se encontraron secuencias en {input_fasta}")

    with open(xml_path, "w") as xml_out, open(summary_path, "w") as summary_out, open(
        text_path, "w"
    ) as text_out:
        summary_out.write("frame\tbest_hit\tevalue\tpct_identity\talign_len\tquery_len\n")
        for rec in sequences:
            print(f"🚀 BLASTP remoto → {rec.id}")

            xml_handle = NCBIWWW.qblast(
                program="blastp",
                database="swissprot",
                sequence=str(rec.seq),
                hitlist_size=hitlist_size,
                format_type="XML",
            )
            xml_text = xml_handle.read()
            xml_handle.close()
            xml_out.write(xml_text)

            blast_record = NCBIXML.read(io.StringIO(xml_text))
            if blast_record.alignments:
                aln = blast_record.alignments[0]
                hsp = aln.hsps[0]
                pct_id = 100.0 * hsp.identities / hsp.align_length
                summary_out.write(
                    f"{rec.id}\t{aln.hit_def}\t{hsp.expect:.2e}\t{pct_id:.2f}\t"
                    f"{hsp.align_length}\t{blast_record.query_length}\n"
                )
            else:
                summary_out.write(
                    f"{rec.id}\t(no hits)\tNA\tNA\t0\t{blast_record.query_length}\n"
                )

            text_handle = NCBIWWW.qblast(
                program="blastp",
                database="swissprot",
                sequence=str(rec.seq),
                hitlist_size=hitlist_size,
                format_type="Text",
            )
            text_out.write(f"# ===== Resultado BLAST para {rec.id} =====\n")
            text_out.write(text_handle.read())
            text_out.write("\n\n")
            text_handle.close()

    print(f"✅ Listo: {xml_path}, {summary_path} y {text_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Ejercicio 2 - BLAST remoto (SwissProt)")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Archivo FASTA generado en el Ejercicio 1")
    parser.add_argument("--out-xml", default=str(DEFAULT_XML), help="Output XML para guardar el BLAST completo")
    parser.add_argument("--summary", default=str(DEFAULT_SUMMARY), help="Resumen TSV con el mejor hit por frame")
    parser.add_argument("--text-output", default=str(DEFAULT_TEXT), help="Archivo de texto estilo blast.out")
    parser.add_argument("--hitlist-size", type=int, default=10, help="Cantidad de hits a recuperar por frame")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_blast_all_frames(
        Path(args.input),
        Path(args.out_xml),
        Path(args.summary),
        Path(args.text_output),
        args.hitlist_size,
    )


if __name__ == "__main__":
    main()
