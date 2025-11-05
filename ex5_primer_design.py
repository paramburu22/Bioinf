#!/usr/bin/env python3
"""
Ejercicio 5 - Diseño Parametrizable de Primers
Diseña primers para el transcript HBB (NM_000518) según parámetros de configuración.
"""

import json
import sys
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import MeltingTemp as mt

# Archivos
GENBANK_FILE = "HBB_NM000518.gb"
CONFIG_FILE = "primer_config.json"
OUTPUT_FASTA = "HBB_primers.fasta"


def load_config(config_file):
    """Carga la configuración desde un archivo JSON."""
    try:
        with open(config_file, "r") as f:
            config = json.load(f)
        return config
    except FileNotFoundError:
        print(f"❌ Archivo de configuración no encontrado: {config_file}")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"❌ Error al leer configuración JSON: {e}")
        sys.exit(1)


def extract_transcript_sequence(genbank_file):
    """Extrae la secuencia de nucleótidos del transcript desde GenBank."""
    try:
        record = SeqIO.read(genbank_file, "genbank")
        seq = record.seq
        print(f"✅ Secuencia extraída: {record.id} ({len(seq)} bp)")
        return str(seq), record.id
    except FileNotFoundError:
        print(f"❌ Archivo GenBank no encontrado: {genbank_file}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error al leer GenBank: {e}")
        sys.exit(1)


def calculate_gc_content(sequence):
    """Calcula el contenido de GC de una secuencia."""
    return gc_fraction(sequence) * 100


def check_terminal_gc(sequence, avoid_positions=2):
    """Verifica si hay G o C en los extremos terminales."""
    if len(sequence) < avoid_positions * 2:
        return False
    
    # Primeros N y últimos N posiciones
    first_n = sequence[:avoid_positions]
    last_n = sequence[-avoid_positions:]
    
    # Verificar si hay G o C en los extremos
    has_gc_start = any(base in 'GC' for base in first_n)
    has_gc_end = any(base in 'GC' for base in last_n)
    
    return not (has_gc_start or has_gc_end)


def calculate_tm(sequence):
    """Calcula la temperatura de melting usando el método de Wallace."""
    try:
        # Usar Tm_Wallace (más simple) o Tm_NN (nearest neighbor)
        # Tm_Wallace: Tm = 64.9 + 41*(GC_count - 16.4) / length
        # Para secuencias cortas, usar Tm_Wallace es más apropiado
        tm = mt.Tm_Wallace(sequence)
        return tm
    except Exception as e:
        # Fallback: cálculo manual simplificado
        gc_count = sequence.count('G') + sequence.count('C')
        length = len(sequence)
        if length == 0:
            return 0
        tm = 64.9 + 41 * (gc_count - 16.4) / length
        return tm


def generate_primer_candidates(sequence, min_length, max_length):
    """Genera candidatos de primers usando ventana deslizante."""
    candidates = []
    
    for length in range(min_length, max_length + 1):
        for i in range(len(sequence) - length + 1):
            primer = sequence[i:i + length]
            candidates.append({
                'sequence': primer,
                'start': i + 1,  # 1-indexed
                'end': i + length,
                'length': length
            })
    
    return candidates


def filter_primers(candidates, config):
    """Filtra primers según los criterios de configuración."""
    filtered = []
    
    min_gc = config.get('min_gc', 50)
    max_gc = config.get('max_gc', 60)
    max_tm = config.get('max_tm', 67)
    avoid_terminal_gc = config.get('avoid_terminal_gc', True)
    terminal_positions = config.get('terminal_positions', 2)
    
    for candidate in candidates:
        seq = candidate['sequence']
        
        # 1. Filtrar por longitud (ya está en el rango)
        if candidate['length'] < config.get('min_length', 18):
            continue
        if candidate['length'] > config.get('max_length', 24):
            continue
        
        # 2. Filtrar por contenido GC
        gc_content = calculate_gc_content(seq)
        if gc_content < min_gc or gc_content > max_gc:
            continue
        
        # 3. Filtrar por GC en extremos terminales
        if avoid_terminal_gc:
            if not check_terminal_gc(seq, terminal_positions):
                continue
        
        # 4. Filtrar por temperatura de melting
        tm = calculate_tm(seq)
        if tm > max_tm:
            continue
        
        # Agregar información adicional
        candidate['gc_content'] = round(gc_content, 2)
        candidate['tm'] = round(tm, 2)
        filtered.append(candidate)
    
    return filtered


def select_top_primers(primers, num_primers):
    """Selecciona los mejores N primers."""
    # Ordenar por temperatura de melting (más cercana a 65°C es mejor)
    # y luego por contenido GC
    target_tm = 65.0
    
    sorted_primers = sorted(
        primers,
        key=lambda x: (abs(x['tm'] - target_tm), abs(x['gc_content'] - 55))
    )
    
    return sorted_primers[:num_primers]


def print_results(primers, transcript_id):
    """Imprime los resultados de los primers diseñados."""
    print("\n" + "=" * 80)
    print("RESULTADOS DEL DISEÑO DE PRIMERS")
    print("=" * 80)
    print(f"Transcript: {transcript_id}")
    print(f"Primers diseñados: {len(primers)}")
    print("=" * 80)
    
    for i, primer in enumerate(primers, 1):
        print(f"\nPrimer {i}:")
        print(f"  Secuencia: {primer['sequence']}")
        print(f"  Posición: {primer['start']}-{primer['end']}")
        print(f"  Longitud: {primer['length']} bp")
        print(f"  Contenido GC: {primer['gc_content']}%")
        print(f"  Temperatura de melting (Tm): {primer['tm']}°C")
    
    print("\n" + "=" * 80)


def save_primers_fasta(primers, transcript_id, output_file):
    """Guarda los primers en formato FASTA."""
    with open(output_file, "w") as f:
        for i, primer in enumerate(primers, 1):
            f.write(f">Primer_{i}_{transcript_id}_pos{primer['start']}-{primer['end']}\n")
            f.write(f"{primer['sequence']}\n")
    
    print(f"✅ Primers guardados en: {output_file}")


def main():
    """Función principal."""
    print("=" * 80)
    print("EJERCICIO 5 - Diseño Parametrizable de Primers")
    print("=" * 80)
    print()
    
    # Cargar configuración
    print(f"📄 Cargando configuración desde {CONFIG_FILE}...")
    config = load_config(CONFIG_FILE)
    print("✅ Configuración cargada")
    print(f"   Longitud: {config.get('min_length')}-{config.get('max_length')} bp")
    print(f"   GC: {config.get('min_gc')}-{config.get('max_gc')}%")
    print(f"   Tm máximo: {config.get('max_tm')}°C")
    print(f"   Número de primers: {config.get('num_primers', 5)}")
    print()
    
    # Extraer secuencia del transcript
    print(f"🔍 Extrayendo secuencia desde {GENBANK_FILE}...")
    sequence, transcript_id = extract_transcript_sequence(GENBANK_FILE)
    print()
    
    # Generar candidatos
    print("🔧 Generando candidatos de primers...")
    candidates = generate_primer_candidates(
        sequence,
        config.get('min_length', 18),
        config.get('max_length', 24)
    )
    print(f"   Candidatos generados: {len(candidates)}")
    
    # Filtrar primers
    print("🔍 Filtrando primers según criterios...")
    filtered = filter_primers(candidates, config)
    print(f"   Primers que cumplen criterios: {len(filtered)}")
    
    if not filtered:
        print("❌ No se encontraron primers que cumplan todos los criterios")
        print("   Considera ajustar los parámetros en el archivo de configuración.")
        sys.exit(1)
    
    # Seleccionar top N primers
    num_primers = config.get('num_primers', 5)
    selected_primers = select_top_primers(filtered, num_primers)
    
    # Mostrar resultados
    print_results(selected_primers, transcript_id)
    
    # Guardar en FASTA
    save_primers_fasta(selected_primers, transcript_id, OUTPUT_FASTA)
    
    print("\n✅ DISEÑO DE PRIMERS COMPLETADO")


if __name__ == "__main__":
    main()

