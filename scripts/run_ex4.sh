#!/bin/bash
# Script para instalar dependencias y ejecutar el Ejercicio 4 - EMBOSS

# No usar set -e para permitir manejo de errores de instalación

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$PROJECT_ROOT"

PYTHON_BIN="${PYTHON_BIN:-python3}"
PIP_BIN="${PIP_BIN:-pip3}"
if [ -n "$VIRTUAL_ENV" ]; then
    if [ -x "$VIRTUAL_ENV/bin/python" ]; then
        PYTHON_BIN="$VIRTUAL_ENV/bin/python"
    fi
    if [ -x "$VIRTUAL_ENV/bin/pip" ]; then
        PIP_BIN="$VIRTUAL_ENV/bin/pip"
    fi
fi

DATA_INTERIM_DIR="data/interim"
RESULTS_EX04="data/results/ex04"
INPUT_FASTA="${DATA_INTERIM_DIR}/hbb_orfs.fasta"
OUTPUT_FILE="${RESULTS_EX04}/hbb_domain_analysis.txt"
EX1_SCRIPT="src/exercises/ex01_generate_orfs.py"
EX4_SCRIPT="src/exercises/ex04_emboss_prosite.py"

echo "================================================================================"
echo "EJERCICIO 4 - Instalación y Ejecución"
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
if ! command_exists "$PYTHON_BIN"; then
    echo -e "${RED}❌ Python no está disponible (${PYTHON_BIN})${NC}"
    echo "Define PYTHON_BIN o activa tu entorno virtual antes de correr este script."
    exit 1
fi
PYTHON_VERSION=$("$PYTHON_BIN" --version)
echo -e "${GREEN}✅ Python encontrado: $PYTHON_VERSION (${PYTHON_BIN})${NC}"
echo ""

# 2. Verificar e instalar Biopython
echo "2. Verificando Biopython..."
if ! "$PYTHON_BIN" -c "import Bio" 2>/dev/null; then
    echo -e "${YELLOW}⚠️  Biopython no está instalado${NC}"
    echo "Instalando Biopython..."
    "$PIP_BIN" install biopython
    echo -e "${GREEN}✅ Biopython instalado${NC}"
else
    echo -e "${GREEN}✅ Biopython ya está instalado${NC}"
fi
echo ""

# 3. Verificar e instalar EMBOSS
echo "3. Verificando EMBOSS..."

# Agregar conda al PATH si está instalado
if [ -d "$HOME/miniconda3" ] && [ ! -z "$(ls -A $HOME/miniconda3/bin/conda 2>/dev/null)" ]; then
    export PATH="$HOME/miniconda3/bin:$PATH"
    # Inicializar conda en el shell actual si es necesario
    eval "$($HOME/miniconda3/bin/conda shell.zsh hook 2>/dev/null || $HOME/miniconda3/bin/conda shell.bash hook 2>/dev/null || true)"
elif [ -d "$HOME/anaconda3" ] && [ ! -z "$(ls -A $HOME/anaconda3/bin/conda 2>/dev/null)" ]; then
    export PATH="$HOME/anaconda3/bin:$PATH"
    eval "$($HOME/anaconda3/bin/conda shell.zsh hook 2>/dev/null || $HOME/anaconda3/bin/conda shell.bash hook 2>/dev/null || true)"
fi

if ! command_exists embossversion; then
    echo -e "${YELLOW}⚠️  EMBOSS no está instalado${NC}"
    
    EMBOSS_INSTALLED=false
    
    # Intentar con conda/bioconda (método recomendado para bioinformática)
    if command_exists conda; then
        echo "Intentando instalar EMBOSS con conda/bioconda..."
        if conda install -c bioconda emboss -y 2>/dev/null; then
            EMBOSS_INSTALLED=true
            echo -e "${GREEN}✅ EMBOSS instalado con conda${NC}"
        fi
    fi
    
    # Si conda falló o no está disponible, intentar otros métodos
    if [ "$EMBOSS_INSTALLED" = false ]; then
        # Detectar sistema operativo
        if [[ "$OSTYPE" == "darwin"* ]]; then
            # macOS - EMBOSS no está en Homebrew estándar
            echo ""
            echo -e "${YELLOW}EMBOSS no está disponible en Homebrew para macOS${NC}"
            echo ""
            echo "Opciones de instalación:"
            echo ""
            echo "1. Instalar con conda/bioconda (RECOMENDADO):"
            echo "   conda install -c bioconda emboss"
            echo ""
            echo "2. Instalar desde la fuente:"
            echo "   Descargar desde: https://emboss.sourceforge.net/download/"
            echo "   O desde GitHub: https://github.com/emboss-dev/emboss"
            echo ""
            echo "3. Usar el script de instalación automática:"
            echo "   ./scripts/install_emboss.sh"
            echo ""
            echo -e "${YELLOW}════════════════════════════════════════════════════════════════════════════${NC}"
            echo -e "${YELLOW}EMBOSS es necesario para ejecutar el Ejercicio 4.${NC}"
            echo ""
            echo "Instalación recomendada:"
            echo "  1. Ejecuta: ./scripts/install_emboss.sh"
            echo "  2. Luego ejecuta este script nuevamente: ./scripts/run_ex4.sh"
            echo ""
            echo -e "${YELLOW}O instala manualmente con conda:${NC}"
            echo "  conda install -c bioconda emboss"
            echo ""
            echo -e "${YELLOW}════════════════════════════════════════════════════════════════════════════${NC}"
            echo ""
            echo -e "${RED}El script no puede continuar sin EMBOSS.${NC}"
            echo "Por favor, instala EMBOSS primero usando una de las opciones arriba."
            exit 1
        elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
            # Linux
            if command_exists apt-get; then
                echo "Instalando EMBOSS con apt-get..."
                sudo apt-get update
                sudo apt-get install -y emboss
                EMBOSS_INSTALLED=true
            elif command_exists yum; then
                echo "Instalando EMBOSS con yum..."
                sudo yum install -y emboss
                EMBOSS_INSTALLED=true
            else
                echo -e "${RED}❌ No se pudo determinar el gestor de paquetes${NC}"
                echo "Instala EMBOSS manualmente desde: https://emboss.sourceforge.net/"
            fi
        fi
        
        if [ "$EMBOSS_INSTALLED" = true ]; then
            echo -e "${GREEN}✅ EMBOSS instalado${NC}"
        fi
    fi
else
    EMBOSS_VERSION=$(embossversion 2>/dev/null | head -1 || echo "instalado")
    echo -e "${GREEN}✅ EMBOSS encontrado: $EMBOSS_VERSION${NC}"
fi
echo ""

# 4. Verificar archivo de entrada del Ejercicio 1
echo "4. Verificando archivo de entrada..."
if [ ! -f "$INPUT_FASTA" ]; then
    echo -e "${YELLOW}⚠️  Archivo ${INPUT_FASTA} no encontrado${NC}"
    echo "Ejecutando Ejercicio 1 para generarlo..."
    if [ -f "$EX1_SCRIPT" ]; then
        "$PYTHON_BIN" "$EX1_SCRIPT"
        if [ ! -f "$INPUT_FASTA" ]; then
            echo -e "${RED}❌ Error al generar ${INPUT_FASTA}${NC}"
            exit 1
        fi
        echo -e "${GREEN}✅ Archivo ${INPUT_FASTA} generado${NC}"
    else
        echo -e "${RED}❌ Script ${EX1_SCRIPT} no encontrado${NC}"
        exit 1
    fi
else
    echo -e "${GREEN}✅ Archivo ${INPUT_FASTA} encontrado${NC}"
fi
echo ""

# 5. Ejecutar el Ejercicio 4
echo "================================================================================"
echo "5. Ejecutando Ejercicio 4..."
echo "================================================================================"
echo ""

if ! "$PYTHON_BIN" "$EX4_SCRIPT"; then
    echo ""
    echo -e "${RED}❌ Error al ejecutar ${EX4_SCRIPT}${NC}"
    exit 1
fi

echo ""
echo "================================================================================"
echo -e "${GREEN}✅ EJERCICIO 4 COMPLETADO${NC}"
echo "================================================================================"
echo "Resultados guardados en: ${OUTPUT_FILE}"

