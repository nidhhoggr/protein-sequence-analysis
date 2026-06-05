#!/bin/bash

# Ebolavirus Analysis Container - Quick Setup Script

set -e

echo "=========================================="
echo "Ebolavirus Analysis Container Setup"
echo "=========================================="
echo ""

# Check if Docker is installed
if ! command -v docker &> /dev/null; then
    echo "❌ Docker is not installed. Please install Docker first."
    exit 1
fi

echo "✓ Docker found: $(docker --version)"

# Check if docker-compose is installed
if ! command -v docker-compose &> /dev/null; then
    echo "⚠ docker-compose not found. Installing with pip..."
    pip install docker-compose
fi

echo "✓ Docker-compose found: $(docker-compose --version)"
echo ""

# Create directories
echo "Creating project directories..."
mkdir -p sequences results scripts logs
echo "✓ Directories created"
echo ""

# Copy Python scripts
echo "Setting up Python scripts..."
if [ -f pipeline.py ]; then
    cp pipeline.py scripts/
    cp utils.py scripts/
    echo "✓ Scripts copied to scripts/"
fi
echo ""

# Build the container
echo "Building Docker image (this may take 10-15 minutes)..."
echo "This will download and install all bioinformatics tools..."
echo ""

docker-compose build

echo ""
echo "=========================================="
echo "✓ Setup Complete!"
echo "=========================================="
echo ""
echo "Next steps:"
echo ""
echo "1. Download ebolavirus sequences:"
echo "   - From NCBI: https://www.ncbi.nlm.nih.gov/genbank/"
echo "   - Or from ViPR: https://www.viprbrc.org/"
echo "   - Save FASTA file to: sequences/"
echo ""
echo "2. Run the analysis pipeline:"
echo "   docker-compose run ebolavirus-analysis \\"
echo "     python scripts/pipeline.py \\"
echo "     -i sequences/your_sequences.fasta \\"
echo "     -o results \\"
echo "     -n 40"
echo ""
echo "3. View results in: results/"
echo ""
echo "For more information, see README.md"
echo ""
