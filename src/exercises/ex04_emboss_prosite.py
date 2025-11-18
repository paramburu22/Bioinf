#!/usr/bin/env python3
"""
Ejercicio 4 - EMBOSS PROSITE Domain Analysis
Analiza dominios de proteínas usando EMBOSS y la base de datos PROSITE.
"""

import os
import subprocess
import sys
import urllib.request
from pathlib import Path

from Bio import SeqIO

# URLs y rutas
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"
DATA_INTERIM = DATA_DIR / "interim"
RESULTS_DIR = DATA_DIR / "results" / "ex04"
TEMP_DIR = RESULTS_DIR / "tmp"
EXTERNAL_DATA = DATA_DIR / "external" / "prosite"
PROSITE_WORKDIR = EXTERNAL_DATA / "emboss_db" / "prosite_db"
PROSITE_DAT_URL = "https://ftp.expasy.org/databases/prosite/prosite.dat"
PROSITE_DOC_URL = "https://ftp.expasy.org/databases/prosite/prosite.doc"
PROSITE_DAT = EXTERNAL_DATA / "prosite.dat"
PROSITE_DOC = EXTERNAL_DATA / "prosite.doc"
INPUT_FASTA = DATA_INTERIM / "hbb_orfs.fasta"
OUTPUT_RESULTS = RESULTS_DIR / "hbb_domain_analysis.txt"


def check_emboss_installed():
    """Verifica si EMBOSS está instalado."""
    try:
        result = subprocess.run(
            ["which", "embossversion"],
            capture_output=True,
            text=True,
            check=True
        )
        print("✅ EMBOSS está instalado")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("❌ EMBOSS no está instalado")
        print("\n📦 Para instalar EMBOSS:")
        print("   macOS: brew install emboss")
        print("   Linux: sudo apt-get install emboss")
        print("   O descargar desde: https://emboss.sourceforge.net/")
        return False


def check_emboss_command(cmd):
    """Verifica si un comando específico de EMBOSS está disponible."""
    try:
        subprocess.run(
            [cmd, "-help"],
            capture_output=True,
            text=True,
            check=True
        )
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def download_prosite():
    """Descarga la base de datos PROSITE si no existe."""
    downloaded = True
    EXTERNAL_DATA.mkdir(parents=True, exist_ok=True)

    if not PROSITE_DAT.exists():
        print("📥 Descargando PROSITE database (prosite.dat)...")
        try:
            urllib.request.urlretrieve(PROSITE_DAT_URL, PROSITE_DAT)
            print("✅ prosite.dat descargado")
        except Exception as exc:  # pylint: disable=broad-except
            print(f"❌ Error al descargar prosite.dat: {exc}")
            downloaded = False

    if not PROSITE_DOC.exists():
        print("📥 Descargando PROSITE documentation (prosite.doc)...")
        try:
            urllib.request.urlretrieve(PROSITE_DOC_URL, PROSITE_DOC)
            print("✅ prosite.doc descargado")
        except Exception as exc:  # pylint: disable=broad-except
            print(f"⚠️  Error al descargar prosite.doc: {exc}")
            print("   prosextract puede fallar sin este archivo")

    return downloaded


def extract_prosite():
    """Prepara la base de datos PROSITE usando EMBOSS prosextract."""
    if not check_emboss_command("prosextract"):
        print("❌ Comando prosextract no disponible")
        return False

    PROSITE_WORKDIR.parent.mkdir(parents=True, exist_ok=True)
    PROSITE_WORKDIR.mkdir(parents=True, exist_ok=True)

    import shutil

    work_dat = PROSITE_WORKDIR / "prosite.dat"
    work_doc = PROSITE_WORKDIR / "prosite.doc"

    if PROSITE_DAT.exists():
        shutil.copy2(PROSITE_DAT, work_dat)
    if PROSITE_DOC.exists():
        shutil.copy2(PROSITE_DOC, work_doc)

    if not work_dat.exists():
        print(f"❌ prosite.dat no encontrado en {work_dat}")
        return False

    if not work_doc.exists():
        print("⚠️  prosite.doc no encontrado - prosextract puede fallar")

    print(f"🔧 Preparando base de datos PROSITE con prosextract en {PROSITE_WORKDIR}...")
    try:
        subprocess.run(
            ["prosextract", "-prositedir", str(PROSITE_WORKDIR)],
            capture_output=True,
            text=True,
            check=True,
        )
        print("✅ Base de datos PROSITE preparada")
        return True
    except subprocess.CalledProcessError as exc:
        print(f"❌ Error al ejecutar prosextract: {exc}")
        print(f"Salida: {exc.stdout}")
        print(f"Error: {exc.stderr}")
        return False


def analyze_domains(fasta_file: Path, output_file: Path):
    """Analiza dominios de proteínas usando EMBOSS patmatmotifs."""
    if not check_emboss_command("patmatmotifs"):
        print("❌ Comando patmatmotifs no disponible")
        return False
    
    print(f"🔍 Analizando dominios en {fasta_file}...")
    
    # Leer todas las secuencias del FASTA
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    
    if not sequences:
        print(f"❌ No se encontraron secuencias en {fasta_file}")
        return False
    
    results = []
    results.append("=" * 80)
    results.append("ANÁLISIS DE DOMINIOS PROSITE - EMBOSS patmatmotifs")
    results.append("=" * 80)
    results.append(f"Archivo de entrada: {fasta_file}")
    results.append(f"Número de secuencias: {len(sequences)}")
    results.append("=" * 80)
    results.append("")
    
    for seq_record in sequences:
        seq_id = seq_record.id.replace('|', '_').replace('/', '_')
        seq = str(seq_record.seq)
        
        results.append(f"\n{'=' * 80}")
        results.append(f"Secuencia: {seq_record.id}")
        results.append(f"Longitud: {len(seq)} aminoácidos")
        results.append(f"{'=' * 80}")
        
        # Crear archivo temporal con una sola secuencia
        TEMP_DIR.mkdir(parents=True, exist_ok=True)
        temp_fasta = TEMP_DIR / f"{seq_id}.fasta"
        patmat_file = TEMP_DIR / f"{seq_id}_patmat.txt"
        
        with open(temp_fasta, "w") as f:
            SeqIO.write(seq_record, f, "fasta")
        
        try:
            # Ejecutar patmatmotifs
            # Nota: patmatmotifs NO acepta -prositedir, usa la base de datos del sistema
            # Si necesitamos usar una base de datos local, debemos configurar el entorno
            # Ejecutar patmatmotifs con opciones para obtener más información
            # -full: documentación completa de los motivos
            # -noprune: incluir patrones simples (modificaciones post-traduccionales)
            cmd = [
                "patmatmotifs",
                "-sequence", str(temp_fasta),
                "-outfile", str(patmat_file),
                "-full",
                "-noprune",
            ]
            
            env = os.environ.copy()
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                env=env
            )
            
            # Leer resultados
            if patmat_file.exists():
                with open(patmat_file, "r") as f:
                    patmat_output = f.read()
                    
                    # Buscar información de dominios en el output
                    lines = patmat_output.split('\n')
                    domain_info = []
                    hit_count = 0
                    current_motif = None
                    in_motif_section = False
                    
                    for i, line in enumerate(lines):
                        line_stripped = line.strip()
                        
                        # Obtener HitCount
                        if 'HitCount:' in line:
                            hit_count = int(line.split('HitCount:')[1].strip())
                        
                        # Buscar secciones de motivos
                        if 'Motif =' in line:
                            current_motif = line.split('Motif =')[1].strip()
                            in_motif_section = True
                            domain_info.append(f"Motivo: {current_motif}")
                        elif 'Length =' in line and in_motif_section:
                            domain_info.append(f"  {line_stripped}")
                        elif ('Start =' in line or 'End =' in line) and in_motif_section:
                            domain_info.append(f"  {line_stripped}")
                        elif line_stripped and not line_stripped.startswith('#') and in_motif_section:
                            # Capturar secuencia del motivo si está presente
                            if len(line_stripped) > 5 and not line_stripped.startswith('-'):
                                if '|' in line_stripped or any(c.isalpha() for c in line_stripped):
                                    domain_info.append(f"  {line_stripped}")
                        elif '---' in line and in_motif_section:
                            in_motif_section = False
                    
                    # Verificar si hay dominios encontrados
                    if hit_count == 0 or not domain_info:
                        results.append("❌ No se encontraron dominios PROSITE")
                    else:
                        results.append(f"✅ Dominios encontrados: {hit_count} motivo(s)")
                        for info in domain_info:
                            results.append(f"  {info}")
            else:
                results.append("⚠️ No se generó archivo de resultados")
                # Mostrar salida estándar si hay
                if result.stdout:
                    results.append(f"Salida: {result.stdout[:500]}")
            
            # Limpiar archivos temporales
            if temp_fasta.exists():
                temp_fasta.unlink()
            if patmat_file.exists():
                patmat_file.unlink()
                
        except subprocess.CalledProcessError as e:
            results.append(f"❌ Error al analizar {seq_record.id}")
            if e.stderr:
                results.append(f"Error: {e.stderr[:200]}")
            if temp_fasta.exists():
                temp_fasta.unlink(missing_ok=True)
            if patmat_file.exists():
                patmat_file.unlink(missing_ok=True)
    
    # Escribir resultados
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as f:
        f.write("\n".join(results))
    
    print(f"✅ Resultados guardados en: {output_file}")
    return True


def main():
    """Función principal."""
    print("=" * 80)
    print("EJERCICIO 4 - EMBOSS PROSITE Domain Analysis")
    print("=" * 80)
    print()
    
    # Verificar EMBOSS
    if not check_emboss_installed():
        sys.exit(1)
    
    # Verificar archivo de entrada
    if not INPUT_FASTA.exists():
        print(f"❌ Archivo de entrada no encontrado: {INPUT_FASTA}")
        print("   Ejecuta primero el Ejercicio 1 para generar este archivo.")
        sys.exit(1)
    
    # Descargar PROSITE
    if not download_prosite():
        sys.exit(1)
    
    # Preparar PROSITE
    if not extract_prosite():
        print("⚠️ Continuando sin prosextract (puede que funcione de todas formas)...")
    
    # Analizar dominios
    if not analyze_domains(INPUT_FASTA, OUTPUT_RESULTS):
        sys.exit(1)
    
    print()
    print("=" * 80)
    print("✅ ANÁLISIS COMPLETADO")
    print("=" * 80)
    print(f"📄 Resultados guardados en: {OUTPUT_RESULTS}")


if __name__ == "__main__":
    main()

