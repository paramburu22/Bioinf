#!/usr/bin/env python3
"""Utilidad para interpretar resultados BLAST tabulados (outfmt 6/7)."""

import argparse
from pathlib import Path
from typing import Dict, Iterable, List


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS_DIR = PROJECT_ROOT / "data" / "results" / "ex02"
DEFAULT_PATTERNS = ("frame_*.out", "blast_frame*.out")


def parse_blast_tab(filename: Path) -> Dict | None:
    """Lee un archivo tabulado (-outfmt 6/7) y devuelve info del mejor hit."""
    with open(filename) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            return {
                "file": filename.name,
                "subject": parts[1],
                "pident": float(parts[2]),
                "length": int(parts[3]),
                "evalue": float(parts[10]),
                "bitscore": float(parts[11]),
                "title": parts[-1].strip(),
            }
    return None


def collect_result_files(base_dir: Path, patterns: Iterable[str]) -> List[Path]:
    files: List[Path] = []
    for pattern in patterns:
        files.extend(sorted(base_dir.glob(pattern)))
    # eliminar duplicados manteniendo orden
    unique_files: List[Path] = []
    seen = set()
    for file in files:
        if file not in seen:
            unique_files.append(file)
            seen.add(file)
    return unique_files


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analiza outputs BLAST tabulados (Ejercicio 2)")
    parser.add_argument("--results-dir", default=str(DEFAULT_RESULTS_DIR), help="Directorio que contiene los .out")
    parser.add_argument(
        "--patterns",
        nargs="+",
        default=list(DEFAULT_PATTERNS),
        help="Patrones glob para detectar outputs (se evalúan dentro de results-dir)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    results_dir = Path(args.results_dir)
    files = collect_result_files(results_dir, args.patterns)
    print(f"🔍 Analizando {len(files)} resultados BLAST dentro de {results_dir}...\n")
    results = []
    for file in files:
        hit = parse_blast_tab(file)
        if hit:
            results.append(hit)

    results.sort(key=lambda x: (x["evalue"], -x["bitscore"]))

    if not results:
        print("⚠️  No se encontraron hits en los archivos especificados.")
        return

    print(f"{'Archivo':<20} {'Top hit':<40} {'%ID':>6} {'E-val':>10}")
    print("-" * 80)
    for r in results:
        print(f"{r['file']:<20} {r['title'][:38]:<40} {r['pident']:>6.1f} {r['evalue']:>10.2e}")

    best = results[0]
    print("\n✅ El marco de lectura correcto es:", best["file"])
    print("   → Mejor hit:", best["title"])
    print("   → % Identidad:", best["pident"], "E-value:", best["evalue"])


if __name__ == "__main__":
    main()
