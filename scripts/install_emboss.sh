#!/bin/bash
# Script para instalar EMBOSS en macOS usando conda/bioconda

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$PROJECT_ROOT"

echo "================================================================================"
echo "Instalación de EMBOSS para macOS"
echo "================================================================================"
echo ""

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Verificar si conda está en el PATH o en ubicaciones comunes
CONDA_PATH=""
if command_exists conda; then
    CONDA_PATH=$(which conda | sed 's|/bin/conda||')
elif [ -f "$HOME/miniconda3/bin/conda" ]; then
    CONDA_PATH="$HOME/miniconda3"
    export PATH="$CONDA_PATH/bin:$PATH"
elif [ -f "$HOME/anaconda3/bin/conda" ]; then
    CONDA_PATH="$HOME/anaconda3"
    export PATH="$CONDA_PATH/bin:$PATH"
fi

# Verificar si conda está disponible ahora
if [ -z "$CONDA_PATH" ] || ! "$CONDA_PATH/bin/conda" --version >/dev/null 2>&1; then
    echo -e "${YELLOW}⚠️  Conda no está instalado o no está en el PATH${NC}"
    echo ""
    
    # Verificar si miniconda3 ya existe
    if [ -d "$HOME/miniconda3" ]; then
        echo -e "${YELLOW}Miniconda ya está instalado en $HOME/miniconda3${NC}"
        echo "Agregando conda al PATH..."
        export PATH="$HOME/miniconda3/bin:$PATH"
        
        # Inicializar conda en el shell actual
        eval "$($HOME/miniconda3/bin/conda shell.zsh hook 2>/dev/null || $HOME/miniconda3/bin/conda shell.bash hook)"
        
        if command_exists conda; then
            echo -e "${GREEN}✅ Conda ahora está disponible${NC}"
            CONDA_PATH="$HOME/miniconda3"
        else
            echo -e "${RED}❌ No se pudo inicializar conda${NC}"
            echo ""
            echo "Por favor, ejecuta manualmente:"
            echo "  export PATH=\"\$HOME/miniconda3/bin:\$PATH\""
            echo "  eval \"\$(\$HOME/miniconda3/bin/conda shell.zsh hook)\""
            echo ""
            echo "O inicializa conda en tu shell:"
            echo "  source ~/.zshrc  # o ~/.bash_profile"
            exit 1
        fi
    else
        echo "Instalando Miniconda..."
        echo ""
        
        # Detectar arquitectura
        ARCH=$(uname -m)
        if [ "$ARCH" = "arm64" ]; then
            MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
        else
            MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
        fi
        
        echo "Descargando Miniconda..."
        curl -O "$MINICONDA_URL"
        
        INSTALLER=$(basename "$MINICONDA_URL")
        echo "Ejecutando instalador..."
        bash "$INSTALLER" -b -p "$HOME/miniconda3"
        
        # Agregar conda al PATH
        export PATH="$HOME/miniconda3/bin:$PATH"
        
        # Inicializar conda
        eval "$($HOME/miniconda3/bin/conda shell.zsh hook 2>/dev/null || $HOME/miniconda3/bin/conda shell.bash hook)"
        
        CONDA_PATH="$HOME/miniconda3"
        echo ""
        echo -e "${GREEN}✅ Miniconda instalado${NC}"
    fi
else
    echo -e "${GREEN}✅ Conda encontrado en: $CONDA_PATH${NC}"
fi

# Asegurar que conda está en el PATH
if [ -n "$CONDA_PATH" ]; then
    export PATH="$CONDA_PATH/bin:$PATH"
    # Inicializar conda en el shell actual
    eval "$($CONDA_PATH/bin/conda shell.zsh hook 2>/dev/null || $CONDA_PATH/bin/conda shell.bash hook 2>/dev/null || true)"
fi

# Verificar que conda funciona
if ! command_exists conda; then
    echo -e "${RED}❌ Conda no está disponible en el PATH${NC}"
    echo "Por favor, inicializa conda manualmente:"
    echo "  export PATH=\"\$HOME/miniconda3/bin:\$PATH\""
    exit 1
fi

echo -e "${GREEN}✅ Conda disponible: $(conda --version)${NC}"
echo ""

# Instalar EMBOSS desde bioconda
echo "Instalando EMBOSS desde bioconda..."
echo ""

# Aceptar términos de servicio de conda si es necesario
echo "Aceptando términos de servicio de conda..."
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 2>/dev/null || true

# Configurar bioconda si no está configurado
if ! conda config --show channels 2>/dev/null | grep -q bioconda; then
    echo "Configurando canal bioconda..."
    conda config --add channels bioconda
    conda config --add channels conda-forge
fi

# Instalar EMBOSS
echo "Instalando EMBOSS..."
conda install -c bioconda emboss -y

echo ""
echo -e "${GREEN}✅ EMBOSS instalado exitosamente${NC}"
echo ""

# Verificar instalación
if command_exists embossversion; then
    EMBOSS_VERSION=$(embossversion 2>/dev/null | head -1 || echo "instalado")
    echo -e "${GREEN}✅ EMBOSS verificado: $EMBOSS_VERSION${NC}"
    echo ""
    echo "Ahora puedes ejecutar el Ejercicio 4:"
    echo "  ./scripts/run_ex4.sh"
else
    echo -e "${YELLOW}⚠️  EMBOSS instalado pero no está en el PATH${NC}"
    echo "Intenta reiniciar tu terminal o ejecutar:"
    echo "  source ~/.bash_profile"
fi

