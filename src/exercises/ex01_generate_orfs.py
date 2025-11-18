#!/usr/bin/env python3
"""Ejercicio 1 - Generación de ORFs desde un GenBank."""

from pathlib import Path

from Bio import SeqIO


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_RAW = PROJECT_ROOT / "data" / "raw"
DATA_INTERIM = PROJECT_ROOT / "data" / "interim"
GENBANK_FILE = DATA_RAW / "hbb_nm000518.gb"
OUTPUT_FASTA = DATA_INTERIM / "hbb_orfs.fasta"


def translate_frames(seq):
    """Genera las seis lecturas posibles (3 directas + 3 reversas complementarias)."""
    frames = []
    for frame in range(3):
        frames.append(seq[frame:].translate(to_stop=False))
    rev = seq.reverse_complement()
    for frame in range(3):
        frames.append(rev[frame:].translate(to_stop=False))
    return frames


def main():
    DATA_INTERIM.mkdir(parents=True, exist_ok=True)
    print(f"🔍 Leyendo {GENBANK_FILE} ...")
    with open(OUTPUT_FASTA, "w") as out_f:
        for record in SeqIO.parse(GENBANK_FILE, "genbank"):
            seq = record.seq
            frames = translate_frames(seq)
            for i, prot in enumerate(frames):
                header = f">{record.id}_frame{i+1}\n"
                out_f.write(header)
                out_f.write(str(prot) + "\n")
    print(f"✅ Archivo FASTA generado: {OUTPUT_FASTA}")


if __name__ == "__main__":
    main()
