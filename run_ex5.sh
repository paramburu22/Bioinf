#!/bin/bash
# Script para instalar dependencias y ejecutar el Ejercicio 5 - Diseño de Primers

set -e  # Salir si hay error

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
if ! command_exists python; then
    echo -e "${RED}❌ python no está instalado${NC}"
    echo "Instala python desde python.org o con: brew install python"
    exit 1
fi
PYTHON_VERSION=$(python --version)
echo -e "${GREEN}✅ Python encontrado: $PYTHON_VERSION${NC}"
echo ""

# 2. Verificar e instalar Biopython
echo "2. Verificando Biopython..."
if ! python -c "import Bio" 2>/dev/null; then
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
if [ ! -f "HBB_NM000518.gb" ]; then
    echo -e "${RED}❌ Archivo HBB_NM000518.gb no encontrado${NC}"
    echo "Este archivo es necesario para extraer la secuencia del transcript."
    exit 1
else
    echo -e "${GREEN}✅ Archivo HBB_NM000518.gb encontrado${NC}"
fi
echo ""

# 4. Verificar archivo de configuración
echo "4. Verificando archivo de configuración..."
if [ ! -f "primer_config.json" ]; then
    echo -e "${YELLOW}⚠️  Archivo primer_config.json no encontrado${NC}"
    echo "Creando archivo de configuración por defecto..."
    cat > primer_config.json << 'EOF'
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
    echo -e "${GREEN}✅ Archivo primer_config.json creado${NC}"
else
    echo -e "${GREEN}✅ Archivo primer_config.json encontrado${NC}"
fi
echo ""

# 5. Ejecutar el Ejercicio 5
echo "================================================================================"
echo "5. Ejecutando Ejercicio 5..."
echo "================================================================================"
echo ""

python ex5_primer_design.py

echo ""
echo "================================================================================"
echo -e "${GREEN}✅ EJERCICIO 5 COMPLETADO${NC}"
echo "================================================================================"
echo "Resultados guardados en: HBB_primers.fasta"

