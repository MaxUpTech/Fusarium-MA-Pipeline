#!/bin/bash
# =============================================================================
# Fusarium MA Pipeline - Dependency Checker and Installer
# =============================================================================
#
# This script checks for all required tools and installs missing ones.
# Run this BEFORE running the pipeline.
#
# Usage: ./check_and_install_dependencies.sh
#
# =============================================================================

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}"
echo "============================================================"
echo "  Fusarium MA Pipeline - Dependency Checker"
echo "============================================================"
echo -e "${NC}"

# Track status
MISSING_TOOLS=()
INSTALLED_TOOLS=()

# Function to check if a command exists
check_tool() {
    local tool=$1
    local install_cmd=$2
    local version_cmd=$3
    
    echo -n "Checking $tool... "
    
    if command -v $tool &> /dev/null; then
        if [ -n "$version_cmd" ]; then
            version=$($version_cmd 2>&1 | head -1)
            echo -e "${GREEN}Found${NC} ($version)"
        else
            echo -e "${GREEN}Found${NC}"
        fi
        INSTALLED_TOOLS+=("$tool")
        return 0
    else
        echo -e "${RED}Not found${NC}"
        MISSING_TOOLS+=("$tool:$install_cmd")
        return 1
    fi
}

# Function to check Python package
check_python_package() {
    local package=$1
    echo -n "Checking Python package $package... "
    
    if python3 -c "import $package" 2>/dev/null; then
        echo -e "${GREEN}Found${NC}"
        return 0
    else
        echo -e "${RED}Not found${NC}"
        return 1
    fi
}

echo -e "${YELLOW}=== Checking Required Tools ===${NC}"
echo ""

# Core bioinformatics tools
echo -e "${BLUE}--- Alignment Tools ---${NC}"
check_tool "bwa-mem2" "conda install -c bioconda bwa-mem2" "bwa-mem2 version"
check_tool "bwa" "conda install -c bioconda bwa" "bwa 2>&1 | grep Version"
check_tool "bowtie2" "conda install -c bioconda bowtie2" "bowtie2 --version"

echo ""
echo -e "${BLUE}--- SAM/BAM Tools ---${NC}"
check_tool "samtools" "conda install -c bioconda samtools" "samtools --version"
check_tool "bcftools" "conda install -c bioconda bcftools" "bcftools --version"

echo ""
echo -e "${BLUE}--- Variant Calling ---${NC}"
check_tool "gatk" "conda install -c bioconda gatk4" "gatk --version"

echo ""
echo -e "${BLUE}--- Quality Control ---${NC}"
check_tool "fastp" "conda install -c bioconda fastp" "fastp --version"
check_tool "fastqc" "conda install -c bioconda fastqc" "fastqc --version"

echo ""
echo -e "${BLUE}--- Annotation ---${NC}"
check_tool "snpEff" "conda install -c bioconda snpeff" "snpEff -version"

echo ""
echo -e "${BLUE}--- Utilities ---${NC}"
check_tool "bedtools" "conda install -c bioconda bedtools" "bedtools --version"
check_tool "vcftools" "conda install -c bioconda vcftools" "vcftools --version"

echo ""
echo -e "${YELLOW}=== Checking Python Environment ===${NC}"
echo ""

# Check Python version
echo -n "Checking Python version... "
PYTHON_VERSION=$(python3 --version 2>&1 | cut -d' ' -f2)
PYTHON_MAJOR=$(echo $PYTHON_VERSION | cut -d'.' -f1)
PYTHON_MINOR=$(echo $PYTHON_VERSION | cut -d'.' -f2)

if [ "$PYTHON_MAJOR" -ge 3 ] && [ "$PYTHON_MINOR" -ge 8 ]; then
    echo -e "${GREEN}$PYTHON_VERSION${NC}"
else
    echo -e "${RED}$PYTHON_VERSION (requires >= 3.8)${NC}"
fi

echo ""
echo -e "${BLUE}--- Python Packages ---${NC}"

PYTHON_PACKAGES=(
    "numpy"
    "scipy"
    "pandas"
    "matplotlib"
    "seaborn"
    "yaml:pyyaml"
    "openpyxl"
    "click"
)

MISSING_PYTHON=()

for pkg in "${PYTHON_PACKAGES[@]}"; do
    # Handle packages with different import names
    if [[ $pkg == *":"* ]]; then
        import_name=$(echo $pkg | cut -d':' -f1)
        pip_name=$(echo $pkg | cut -d':' -f2)
    else
        import_name=$pkg
        pip_name=$pkg
    fi
    
    if ! check_python_package "$import_name"; then
        MISSING_PYTHON+=("$pip_name")
    fi
done

# Summary
echo ""
echo -e "${YELLOW}=== Summary ===${NC}"
echo ""

if [ ${#MISSING_TOOLS[@]} -eq 0 ] && [ ${#MISSING_PYTHON[@]} -eq 0 ]; then
    echo -e "${GREEN}All dependencies are installed!${NC}"
    echo ""
    echo "You can now run the pipeline:"
    echo "  python3 scripts/run_pipeline.py -c config/your_config.yaml"
    exit 0
fi

# Report missing tools
if [ ${#MISSING_TOOLS[@]} -gt 0 ]; then
    echo -e "${RED}Missing bioinformatics tools:${NC}"
    for item in "${MISSING_TOOLS[@]}"; do
        tool=$(echo $item | cut -d':' -f1)
        cmd=$(echo $item | cut -d':' -f2)
        echo "  - $tool"
    done
    echo ""
fi

if [ ${#MISSING_PYTHON[@]} -gt 0 ]; then
    echo -e "${RED}Missing Python packages:${NC}"
    for pkg in "${MISSING_PYTHON[@]}"; do
        echo "  - $pkg"
    done
    echo ""
fi

# Offer to install
echo -e "${YELLOW}Would you like to install missing dependencies? (y/n)${NC}"
read -r response

if [[ "$response" =~ ^[Yy]$ ]]; then
    echo ""
    echo -e "${BLUE}Installing dependencies...${NC}"
    echo ""
    
    # Install Python packages
    if [ ${#MISSING_PYTHON[@]} -gt 0 ]; then
        echo "Installing Python packages..."
        pip3 install ${MISSING_PYTHON[@]}
    fi
    
    # Install bioinformatics tools via conda
    if [ ${#MISSING_TOOLS[@]} -gt 0 ]; then
        echo ""
        echo "Installing bioinformatics tools via conda..."
        echo "(Make sure conda is activated)"
        
        CONDA_PACKAGES=""
        for item in "${MISSING_TOOLS[@]}"; do
            tool=$(echo $item | cut -d':' -f1)
            case $tool in
                bwa-mem2) CONDA_PACKAGES+=" bwa-mem2" ;;
                bwa) CONDA_PACKAGES+=" bwa" ;;
                bowtie2) CONDA_PACKAGES+=" bowtie2" ;;
                samtools) CONDA_PACKAGES+=" samtools" ;;
                bcftools) CONDA_PACKAGES+=" bcftools" ;;
                gatk) CONDA_PACKAGES+=" gatk4" ;;
                fastp) CONDA_PACKAGES+=" fastp" ;;
                fastqc) CONDA_PACKAGES+=" fastqc" ;;
                snpEff) CONDA_PACKAGES+=" snpeff" ;;
                bedtools) CONDA_PACKAGES+=" bedtools" ;;
                vcftools) CONDA_PACKAGES+=" vcftools" ;;
            esac
        done
        
        if [ -n "$CONDA_PACKAGES" ]; then
            echo "Running: conda install -c bioconda -c conda-forge $CONDA_PACKAGES"
            conda install -c bioconda -c conda-forge -y $CONDA_PACKAGES
        fi
    fi
    
    echo ""
    echo -e "${GREEN}Installation complete!${NC}"
    echo "Please run this script again to verify all dependencies."
else
    echo ""
    echo "To install manually:"
    echo ""
    
    if [ ${#MISSING_PYTHON[@]} -gt 0 ]; then
        echo "  pip3 install ${MISSING_PYTHON[@]}"
    fi
    
    if [ ${#MISSING_TOOLS[@]} -gt 0 ]; then
        echo "  conda install -c bioconda -c conda-forge bwa-mem2 samtools bcftools gatk4 fastp fastqc snpeff bedtools"
    fi
fi
