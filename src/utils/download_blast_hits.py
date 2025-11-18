#!/usr/bin/env python3
"""Descarga los top-N hits de un resultado BLAST (XML) y los guarda en FASTA."""

import argparse
import os
from pathlib import Path

from Bio import Entrez
from Bio.Blast import NCBIXML


PROJECT_ROOT = Path(__file__).resolve().parents[2]
EX02_RESULTS = PROJECT_ROOT / "data" / "results" / "ex02"
EX03_RESULTS = PROJECT_ROOT / "data" / "results" / "ex03"


def download_hits(blast_xml: Path, output_fasta: Path, email: str, top_n: int, api_key: str | None = None) -> None:
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    print(f"🔍 Leyendo resultados de BLAST → {blast_xml}")
    with open(blast_xml) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        try:
            record = next(blast_records)
        except StopIteration as exc:
            raise SystemExit("El archivo BLAST está vacío.") from exc

    hits = [alignment.accession for alignment in record.alignments[:top_n]]
    if not hits:
        raise SystemExit("No se encontraron hits en el archivo BLAST.")

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    print("\n⬇️  Descargando secuencias desde NCBI...")
    with open(output_fasta, "w") as out_f:
        for acc in hits:
            print(f"   → {acc}")
            with Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text") as handle:
                out_f.write(handle.read())
    print(f"\n✅ Archivo final creado: {output_fasta}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Descarga los top hits desde un BLAST XML.")
    parser.add_argument("--blast-xml", default=str(EX02_RESULTS / "blast_hbb_remote.xml"), help="Resultado BLAST en formato XML")
    parser.add_argument("--output-fasta", default=str(EX03_RESULTS / "msa_input.fasta"), help="Archivo FASTA para los hits descargados")
    parser.add_argument("--top", type=int, default=10, help="Número de hits a descargar")
    parser.add_argument("--email", default=os.environ.get("ENTREZ_EMAIL"), help="Email requerido por NCBI (o variable ENTREZ_EMAIL)")
    parser.add_argument("--api-key", default=os.environ.get("ENTREZ_API_KEY"), help="API key opcional de NCBI")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.email:
        raise SystemExit("Es obligatorio especificar un email (--email o variable ENTREZ_EMAIL).")
    download_hits(Path(args.blast_xml), Path(args.output_fasta), args.email, args.top, args.api_key)


if __name__ == "__main__":
    main()
