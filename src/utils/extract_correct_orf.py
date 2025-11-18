#!/usr/bin/env python3
"""Extrae el CDS principal (feature CDS) desde un GenBank y lo traduce a proteína."""

import argparse
from pathlib import Path

from Bio import SeqIO


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_RAW = PROJECT_ROOT / "data" / "raw"
DATA_INTERIM = PROJECT_ROOT / "data" / "interim"


def extract_correct_orf(genbank_file: Path, output_fasta: Path, header: str = "HBB_protein_correct_ORF") -> None:
    record = SeqIO.read(genbank_file, "genbank")
    cds_feats = [f for f in record.features if f.type == "CDS"]
    if not cds_feats:
        raise SystemExit("No se encontró feature CDS en el GenBank.")

    cds = cds_feats[0]
    cds_seq = cds.extract(record.seq)
    prot = cds_seq.translate(to_stop=True)

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(output_fasta, "w") as f:
        f.write(f">{header}\n{prot}\n")
    print(f"Listo → {output_fasta} len={len(prot)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extrae y traduce el CDS principal de un GenBank.")
    parser.add_argument("--genbank", default=str(DATA_RAW / "hbb_nm000518.gb"), help="Archivo GenBank de entrada")
    parser.add_argument(
        "--output",
        default=str(DATA_INTERIM / "hbb_correct_orf.fasta"),
        help="Archivo FASTA de salida con la proteína traducida",
    )
    parser.add_argument("--header", default="HBB_protein_correct_ORF", help="Encabezado FASTA para la proteína")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    extract_correct_orf(Path(args.genbank), Path(args.output), args.header)


if __name__ == "__main__":
    main()
