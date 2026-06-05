# Protein Sequence Diversity Analysis

A complete containerized bioinformatics pipeline for downloading protein sequences from UniProt (multi species support), removing redundancy, aligning, building phylogenetic trees, and selecting diverse representatives across all species to avoid sampling bias.

## Features

✅ **UniProt integration** - Direct downloads of curated sequences  
✅ **Redundancy removal** - CD-HIT clustering at 95% identity threshold  
✅ **Multiple alignment** - MAFFT optimized for protein sequences  
✅ **Phylogenetic inference** - IQ-TREE with 1000 bootstrap replicates  
✅ **Diversity-aware selection** - Phylogenetic tree-guided sampling  
✅ **No species bias** - Balanced representation across all 6 ebolavirus species  
✅ **Fully configurable** - JSON config files for species, proteins, and parameters  
✅ **Complete logging** - Full pipeline logs and statistics  

## System Requirements

- **Docker** (with docker compose)
- **4GB RAM** minimum (8GB+ recommended)
- **20GB disk space** (for sequences + intermediate files)
- **Internet connection** (for UniProt downloads)

## Installation

### 1. Prerequisites

Make sure Docker is installed:
```bash
docker --version
docker compose version
```

### 2. Setup Project

```bash
git clone git@github.com:nidhhoggr/protein-sequence-analysis.git ebolavirus-analysis
```

### 3. Build Container

```bash
docker compose build
```

This takes 10-15 minutes the first time (installs all bioinformatics tools).

## Quick Start

### Step 1: Download Glycoprotein Sequences

Download all 6 species from UniProt:

```bash
docker compose run ebolavirus-analysis \
  python scripts/download_gp_uniprot.py \
  -o sequences
```

**Output:** 
- Individual files for each species:
  - `sequences/ebolavirus_gp_zaire.fasta`
  - `sequences/ebolavirus_gp_sudan.fasta`
  - `sequences/ebolavirus_gp_bundibugyo.fasta`
  - `sequences/ebolavirus_gp_tai_forest.fasta`
  - `sequences/ebolavirus_gp_reston.fasta`
  - `sequences/ebolavirus_gp_bombali.fasta`
- Combined file: `sequences/ebolavirus_gp_combined.fasta`
- Total: ~550 sequences across all species

### Step 2: Analyze Species Coverage

Check that you have sequences from all 6 species:

```bash
docker compose run ebolavirus-analysis \
  python scripts/analyze_subtypes.py \
  sequences/ebolavirus_gp_combined.fasta
```

This shows counts for each species and warns if any are missing.

### Step 3: Run Full Analysis Pipeline

This performs clustering, alignment, tree building, and diversity selection:

```bash
docker compose run ebolavirus-analysis \
  python scripts/pipeline.py \
  -i sequences/ebolavirus_gp_combined.fasta \
  -o results \
  -n 40
```

**Pipeline steps:**
1. **CD-HIT clustering** (5 min) - Removes redundant sequences at 95% identity
2. **MAFFT alignment** (2-5 min) - Aligns protein sequences
3. **IQ-TREE inference** (10-30 min) - Builds phylogenetic tree with bootstraps
4. **Phylogenetic selection** (<1 min) - Picks 40 diverse sequences

**Total runtime:** ~30-45 minutes

### Step 4: View Results

```bash
# List output files
ls -lh results/

# Count final sequences
grep "^>" results/ebolavirus_diverse_dataset.fasta | wc -l

# View sequence IDs
grep "^>" results/ebolavirus_diverse_dataset.fasta | head -20
```

## Output Files

After running the pipeline, `results/` contains:

| File | Description |
|------|-------------|
| `ebolavirus_diverse_dataset.fasta` | **Your final result**: 40 diverse sequences |
| `sequences_nonredundant.fasta` | After CD-HIT clustering |
| `sequences_aligned.fasta` | MAFFT multiple sequence alignment |
| `tree.treefile` | Newick format phylogenetic tree |
| `tree.iqtree` | Detailed IQ-TREE output & statistics |
| `pipeline.log` | Full execution log |

## Configurable Settings

All three main scripts support custom configuration via JSON config files.

### Configuration File Format

Create a config file (e.g., `config/config.json`):

```json
{
  "species": {
    "Zaire": "organism_name:Zaire ebolavirus",
    "Sudan": "organism_name:Sudan ebolavirus",
    "Bundibugyo": "organism_name:Bundibugyo ebolavirus",
    "Taï Forest": "organism_name:Taï Forest ebolavirus",
    "Reston": "organism_name:Reston ebolavirus",
    "Bombali": "organism_name:Bombali virus"
  },
  "protein_filter": "protein_name:glycoprotein"
}
```

**For download_gp_uniprot.py and pipeline.py:**
- `species`: UniProt organism_name queries for each species
- `protein_filter`: UniProt protein_name filter (glycoprotein, nucleoprotein, polymerase, etc.)

**For analyze_subtypes.py:**
```json
{
  "Zaire": ["Zaire ebolavirus", "EBOV"],
  "Sudan": ["Sudan ebolavirus", "SUDV"],
  "Bundibugyo": ["Bundibugyo", "BDBV"],
  "Taï Forest": ["Taï Forest", "TAFV"],
  "Reston": ["Reston", "RESTV"],
  "Bombali": ["Bombali", "BOMV"]
}
```

### Download Proteins

Then:
```bash
docker compose run ebolavirus-analysis \
  python scripts/download_gp_uniprot.py \
  -o sequences_np \
  -c config/config.json
```

### Use Custom Config in Pipeline

**Download script:**
```bash
docker compose run ebolavirus-analysis \
  python scripts/download_gp_uniprot.py \
  -o sequences \
  -c config/config.json
```

**Analysis script:**
```bash
docker compose run ebolavirus-analysis \
  python scripts/analyze_subtypes.py \
  sequences/ebolavirus_gp_combined.fasta \
  -c config/analyze_config.json
```

**Pipeline script:**
```bash
docker compose run ebolavirus-analysis \
  python scripts/pipeline.py \
  -i sequences/ebolavirus_gp_combined.fasta \
  -o results \
  -n 40 \
  -c config/config.json
```

## Visualizing Results

### View the Phylogenetic Tree

**Online (easiest):**
1. Go to [iTOL](https://itol.embl.de/)
2. Upload `results/tree.treefile`
3. Customize colors, labels, scale

**Locally:**
```bash
# Install FigTree
brew install figtree  # macOS
# Then open the tree file
figtree results/tree.treefile
```

### View the Alignment

Download and use [AliView](http://www.ormbunkar.se/aliview/):
```bash
aliview results/sequences_aligned.fasta
```

## Advanced Usage

### Change Number of Selected Sequences

```bash
docker compose run ebolavirus-analysis \
  python scripts/pipeline.py \
  -i sequences/ebolavirus_gp_combined.fasta \
  -o results \
  -n 50  # Select 50 instead of 40
```

### Use Faster Tree Building

For large datasets, use FastTree instead of IQ-TREE:

```bash
docker compose run ebolavirus-analysis \
  python scripts/pipeline.py \
  -i sequences/ebolavirus_gp_combined.fasta \
  -o results \
  -n 40 \
  --use-fasttree
```

### Interactive Session

For manual analysis or debugging:

```bash
docker compose run ebolavirus-analysis bash
```

Now you have access to all tools:
```bash
# Inside container:
cd /data
mafft --version
iqtree --version
cd-hit --version
python3
```

### Run Individual Steps

```bash
# Just clustering
docker compose run ebolavirus-analysis bash -c \
  "cd-hit -i sequences/ebolavirus_gp_combined.fasta \
   -o results/nr.fasta -c 0.95 -n 5"

# Just alignment
docker compose run ebolavirus-analysis bash -c \
  "mafft --auto results/nr.fasta > results/aligned.fasta"

# Just tree building
docker compose run ebolavirus-analysis bash -c \
  "iqtree -s results/aligned.fasta -st AA -m LG+G -bb 1000 -nt AUTO"
```

## Pipeline Details

### What Each Species Contributes

From UniProt downloads:
- **Zaire ebolavirus**: 487 sequences (most common)
- **Sudan ebolavirus**: 24 sequences
- **Bundibugyo ebolavirus**: 3 sequences
- **Taï Forest ebolavirus**: 11 sequences
- **Reston ebolavirus**: 25 sequences
- **Bombali virus**: 8 sequences

### Diversity Selection Strategy

The pipeline uses phylogenetic guidance to avoid bias:

1. **Clustering** removes near-identical duplicates (keeps structural variants)
2. **Tree building** reveals evolutionary relationships between sequences
3. **Clade-aware sampling** selects representatives from each major branch
4. **Result**: Balanced dataset across all species, no oversampling of similar strains

### Quality Metrics

The pipeline reports:
- Total sequences processed
- Sequences after redundancy removal
- Alignment length and gaps
- Tree depth and clade structure
- Final diversity statistics

## Troubleshooting

### "docker compose: command not found"

Docker Compose should come with Docker Desktop (macOS/Windows). On Linux:
```bash
sudo apt-get install docker-compose-plugin
```

Or use older syntax:
```bash
docker-compose run ...  # Old syntax
docker compose run ...  # New syntax (preferred)
```

### Out of Memory Errors

Increase Docker's memory allocation:

**Docker Desktop GUI:**
- Settings → Resources → Memory → Increase to 8GB+

**Or edit docker-compose.yml:**
```yaml
services:
  ebolavirus-analysis:
    mem_limit: 8g
```

### IQ-TREE Takes Too Long

Use FastTree instead (faster, slightly less accurate):
```bash
docker compose run ebolavirus-analysis \
  python scripts/pipeline.py \
  -i sequences/ebolavirus_gp_combined.fasta \
  -o results --use-fasttree
```

### Sequences Not Found

Make sure the file path is correct:
```bash
docker compose run ebolavirus-analysis bash -c \
  "ls -lh sequences/"
```

Files should be in `sequences/` directory directly.

### Pipeline Crashes Midway

Check the log file:
```bash
docker compose run ebolavirus-analysis bash -c \
  "tail -100 results/pipeline.log"
```

Common issues:
- Out of disk space: `df -h`
- Out of memory: Increase Docker memory
- Network issues: Check internet connection

## Clean Up

### Remove Orphan Containers

```bash
docker compose down --remove-orphans
```

### Delete Results and Start Over

```bash
rm -rf results/*
docker compose run ebolavirus-analysis \
  python scripts/pipeline.py \
  -i sequences/ebolavirus_gp_combined.fasta \
  -o results -n 40
```

### Clean Up Everything

```bash
docker compose down -v
rm -rf results sequences
```

## Next Steps with Your Data

Once you have your 40 diverse sequences:

1. **Structural analysis** - Use ColabFold for 3D protein structure prediction
2. **Sequence evolution** - Analyze conservation with ConSurf
3. **Antibody targeting** - Identify epitopes for vaccine design
4. **Molecular dynamics** - Simulate protein behavior with GROMACS
5. **Comparative genomics** - Study functional domains across species

## Citations

If you use this pipeline, please cite:

- **MAFFT**: Katoh K, et al. (2005) MAFFT: a novel method for rapid multiple sequence alignment
- **IQ-TREE**: Nguyen LT, Schmidt HA, et al. (2015) IQ-TREE: A fast and effective stochastic algorithm
- **CD-HIT**: Li W, Godzik A. (2006) Cd-hit: a fast program for clustering and comparing large sets
- **Biopython**: Cock PJ, et al. (2009) Biopython: freely available Python tools for computational molecular biology

## Support

For issues:

1. Check `results/pipeline.log` for error details
2. Review this README's troubleshooting section
3. Verify Docker installation: `docker compose version`
4. Check available disk/memory: `df -h` and `free -h`
5. Test UniProt connectivity: `curl https://rest.uniprot.org/uniprotkb/search`

## License

This pipeline is provided for research purposes. See individual tool licenses for terms.

---

**Happy analyzing!** 🧬
