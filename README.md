# Bioinf - Trabajo Práctico

Repositorio para ejercicios de Bioinformática.

## Estructura del Proyecto

- `ex1.py` - Ejercicio 1: Generación de ORFs desde GenBank
- `ex2.py` - Ejercicio 2: Análisis BLAST
- `ex4_emboss_analysis.py` - Ejercicio 4: Análisis de dominios PROSITE con EMBOSS
- `ex5_primer_design.py` - Ejercicio 5: Diseño parametrizable de primers
- `run_ex4.sh` - Script de instalación y ejecución del Ejercicio 4
- `run_ex5.sh` - Script de instalación y ejecución del Ejercicio 5
- `install_emboss.sh` - Script auxiliar para instalar EMBOSS con conda/bioconda
- `primer_config.json` - Archivo de configuración para diseño de primers

## Requisitos Previos

### Para el Ejercicio 4:
- Python 3.x
- Biopython
- EMBOSS

### Para el Ejercicio 5:
- Python 3.x
- Biopython

## Instalación y Ejecución

### Ejercicio 4 - Análisis EMBOSS PROSITE

**IMPORTANTE:** EMBOSS no está disponible en Homebrew para macOS. Se recomienda usar conda/bioconda.

#### Opción 1: Instalación automática con conda (RECOMENDADO)

```bash
# Si tienes conda instalado:
conda install -c bioconda emboss

# Si no tienes conda, usa el script de instalación:
./install_emboss.sh
```

Luego ejecuta el ejercicio:

```bash
./run_ex4.sh
```

#### Opción 2: Instalación manual

```bash
# Instalar dependencias Python
pip install biopython

# Instalar EMBOSS (requiere conda):
conda install -c bioconda emboss

# O desde la fuente:
# Descargar desde: https://emboss.sourceforge.net/download/
# O desde GitHub: https://github.com/emboss-dev/emboss

# Ejecutar ejercicio
python ex4_emboss_analysis.py
```

**Nota:** En Linux, EMBOSS puede instalarse con:
```bash
sudo apt-get install emboss  # Debian/Ubuntu
# o
sudo yum install emboss      # RHEL/CentOS
```

**Output:** `HBB_domain_analysis.txt` - Archivo con resultados del análisis de dominios PROSITE.

### Ejercicio 5 - Diseño de Primers

El script `run_ex5.sh` instala automáticamente todas las dependencias y ejecuta el ejercicio:

```bash
./run_ex5.sh
```

O manualmente:

```bash
# Instalar dependencias
pip install biopython

# Ejecutar ejercicio
python ex5_primer_design.py
```

**Output:** 
- `HBB_primers.fasta` - Archivo FASTA con los 5 primers diseñados
- Resultados en consola con detalles de cada primer

### Configuración del Ejercicio 5

El archivo `primer_config.json` contiene los parámetros de diseño:

```json
{
  "min_length": 18,
  "max_length": 24,
  "min_gc": 50,
  "max_gc": 60,
  "max_tm": 67,
  "num_primers": 5,
  "avoid_terminal_gc": true,
  "terminal_positions": 2
}
```

## Archivos de Entrada

- `HBB_NM000518.gb` - Archivo GenBank con el transcript HBB
- `HBB_ORFs.fasta` - Output del Ejercicio 1 (generado automáticamente si no existe)

## Dependencias

Instalar todas las dependencias:

```bash
pip3 install -r requirements.txt
```

## Notas

- El Ejercicio 4 requiere que el Ejercicio 1 haya sido ejecutado previamente para generar `HBB_ORFs.fasta`
- El Ejercicio 5 requiere el archivo `HBB_NM000518.gb` para extraer la secuencia del transcript
- Los scripts de instalación detectan automáticamente el sistema operativo (macOS/Linux)
