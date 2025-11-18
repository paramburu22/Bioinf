#!/bin/bash
# Script para instalar dependencias y ejecutar el Ejercicio 5 - Diseño de Primers

set -e  # Salir si hay error

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$PROJECT_ROOT"

GENBANK_FILE="data/raw/hbb_nm000518.gb"
CONFIG_FILE="config/primer_design.json"
RESULTS_FASTA="data/results/ex05/hbb_primers.fasta"
EX5_SCRIPT="src/exercises/ex05_primer_design.py"

echo "================================================================================"
echo "EJERCICIO 5 - Instalación y Ejecución"
echo "================================================================================"
echo ""

# Colores para output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Función para verificar si un comando existe
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# 1. Verificar Python
echo "1. Verificando Python..."
if ! command_exists python3; then
    echo -e "${RED}❌ python3 no está instalado${NC}"
    echo "Instala python3 desde python.org o con: brew install python"
    exit 1
fi
PYTHON_VERSION=$(python3 --version)
echo -e "${GREEN}✅ Python encontrado: $PYTHON_VERSION${NC}"
echo ""

# 2. Verificar e instalar Biopython
echo "2. Verificando Biopython..."
if ! python3 -c "import Bio" 2>/dev/null; then
    echo -e "${YELLOW}⚠️  Biopython no está instalado${NC}"
    echo "Instalando Biopython..."
    pip3 install biopython
    echo -e "${GREEN}✅ Biopython instalado${NC}"
else
    echo -e "${GREEN}✅ Biopython ya está instalado${NC}"
fi
echo ""

# 3. Verificar archivo GenBank
echo "3. Verificando archivo GenBank..."
if [ ! -f "$GENBANK_FILE" ]; then
    echo -e "${RED}❌ Archivo ${GENBANK_FILE} no encontrado${NC}"
    echo "Este archivo es necesario para extraer la secuencia del transcript."
    exit 1
else
    echo -e "${GREEN}✅ Archivo ${GENBANK_FILE} encontrado${NC}"
fi
echo ""

# 4. Verificar archivo de configuración
echo "4. Verificando archivo de configuración..."
if [ ! -f "$CONFIG_FILE" ]; then
    echo -e "${YELLOW}⚠️  Archivo ${CONFIG_FILE} no encontrado${NC}"
    echo "Creando archivo de configuración por defecto..."
    mkdir -p "$(dirname "$CONFIG_FILE")"
    cat > "$CONFIG_FILE" << 'EOF'
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
EOF
    echo -e "${GREEN}✅ Archivo ${CONFIG_FILE} creado${NC}"
else
    echo -e "${GREEN}✅ Archivo ${CONFIG_FILE} encontrado${NC}"
fi
echo ""

# 5. Ejecutar el Ejercicio 5
echo "================================================================================"
echo "5. Ejecutando Ejercicio 5..."
echo "================================================================================"
echo ""

python3 "$EX5_SCRIPT"

echo ""
echo "================================================================================"
echo -e "${GREEN}✅ EJERCICIO 5 COMPLETADO${NC}"
echo "================================================================================"
echo "Resultados guardados en: ${RESULTS_FASTA}"

