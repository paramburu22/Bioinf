# Ex1.py
from Bio import SeqIO

def translate_frames(seq):
    """Genera las seis lecturas posibles (3 directas + 3 reversas complementarias)."""
    frames = []
    # 3 marcos de lectura directos
    for frame in range(3):
        frames.append(seq[frame:].translate(to_stop=False))
    # 3 marcos de lectura reversos
    rev = seq.reverse_complement()
    for frame in range(3):
        frames.append(rev[frame:].translate(to_stop=False))
    return frames

def main():
    input_file = "HBB_NM000518.gb"
    output_file = "HBB_ORFs.fasta"

    print(f"🔍 Leyendo {input_file} ...")
    with open(output_file, "w") as out_f:
        for record in SeqIO.parse(input_file, "genbank"):
            seq = record.seq
            frames = translate_frames(seq)
            for i, prot in enumerate(frames):
                header = f">{record.id}_frame{i+1}\n"
                out_f.write(header)
                out_f.write(str(prot) + "\n")

    print(f"✅ Archivo FASTA generado: {output_file}")

if __name__ == "__main__":
    main()
