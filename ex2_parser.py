# Ex2_parser.py
import glob

def parse_blast_tab(filename):
    """Lee un archivo de salida BLAST tabulado (-outfmt 6/7) y devuelve info del mejor hit."""
    with open(filename) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            return {
                "file": filename,
                "subject": parts[1],
                "pident": float(parts[2]),
                "length": int(parts[3]),
                "evalue": float(parts[10]),
                "bitscore": float(parts[11]),
                "title": parts[-1].strip(),
            }
    return None

def main():
    files = sorted(glob.glob("frame_*.out")) or sorted(glob.glob("blast_frame*.out"))
    print(f"🔍 Analizando {len(files)} resultados BLAST...\n")
    results = []
    for f in files:
        hit = parse_blast_tab(f)
        if hit:
            results.append(hit)

    # Ordenar por E-value y bitscore
    results.sort(key=lambda x: (x["evalue"], -x["bitscore"]))

    print(f"{'Archivo':<20} {'Top hit':<40} {'%ID':>6} {'E-val':>10}")
    print("-" * 80)
    for r in results:
        print(f"{r['file']:<20} {r['title'][:38]:<40} {r['pident']:>6.1f} {r['evalue']:>10.2e}")

    # Mostrar el mejor frame (menor E-value)
    if results:
        best = results[0]
        print("\n✅ El marco de lectura correcto es:", best["file"])
        print("   → Mejor hit:", best["title"])
        print("   → % Identidad:", best["pident"], "E-value:", best["evalue"])

if __name__ == "__main__":
    main()
