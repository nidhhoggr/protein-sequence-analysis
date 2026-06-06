# Protein Sequence Diversity Analysis

A complete containerized bioinformatics pipeline for downloading protein sequences from UniProt, removing redundancy, aligning, building phylogenetic trees, and selecting diverse representatives while avoiding sampling bias.

**Works for any organism** - not just ebolavirus. Choose between phylogenetic or species-balanced diversity selection.

## Features

✅ **UniProt integration** - Direct downloads of curated sequences  
✅ **Redundancy removal** - CD-HIT clustering (configurable threshold)  
✅ **Multiple alignment** - MAFFT optimized for protein/nucleotide sequences  
✅ **Phylogenetic inference** - IQ-TREE with 1000 bootstrap replicates or FastTree  
✅ **Two diversity strategies**:
  - Tree-guided selection (maximum evolutionary coverage)
  - Species-balanced selection (equal per species)  
✅ **No sampling bias** - Balanced representation across all species  
✅ **Fully configurable** - JSON config files and command-line arguments  
✅ **Generic design** - Works for any organism, not just viruses  
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
git clone git@github.com:nidhhoggr/protein-sequence-analysis.git
cd protein-sequence-analysis
```

### 3. Build Container

```bash
docker compose build
```

Takes 10-15 minutes the first time (installs all bioinformatics tools).

## Quick Start

### Step 1: Download Sequences from UniProt

For ebolavirus example:
```bash
docker compose run ebolavirus-analysis \
  python scripts/download_gp_uniprot.py \
  -o sequences
```

For custom organism, edit `config/config.json` with your species and run above.

**Output:**
- Individual files per species in `sequences/`
- Combined file: `sequences/ebolavirus_gp_combined.fasta`
- ~460 total sequences (if using envelope glycoprotein filter)

### Step 2: Choose Your Diversity Strategy

#### Option A: Species-Balanced Selection (Recommended for imbalanced data)

Use this when you have many sequences from one species and few from others:

```bash
docker compose run ebolavirus-analysis \
  python scripts/balanced_diversity_pipeline.py \
  -i sequences/ebolavirus_gp_combined.fasta \
  --target-per-species 5
```

**What it does:**
1. Selects equal number of sequences from each species (balanced)
2. CD-HIT clustering to remove near-duplicates
3. MAFFT alignment
4. IQ-TREE tree building

**Best for:** Imbalanced datasets (e.g., 400 Zaire vs 3 Bundibugyo)

#### Option B: Phylogenetic Tree-Guided Selection (Maximum diversity)

Use this when you want maximum evolutionary coverage:

```bash
docker compose run ebolavirus-analysis \
  python scripts/generic_pipeline.py \
  -i sequences/ebolavirus_gp_combined.fasta \
  -n 40
```

**What it does:**
1. CD-HIT clustering to remove redundancy
2. MAFFT alignment
3. IQ-TREE tree building
4. Selects 40 sequences spread across phylogenetic tree

**Best for:** Large datasets where you want diversity across all clades

### Step 3: View Results

```bash
# List output files
ls -lh results/

# For balanced approach
grep "^>" results/sequences_balanced.fasta | wc -l
grep "^>" results/sequences_clustered.fasta | wc -l
grep "^>" results/sequences_aligned.fasta | wc -l

# For tree-guided approach  
grep "^>" results/diverse_dataset.fasta | wc -l
```

### Step 4: Visualize Phylogenetic Tree

**Online (easiest):**
1. Go to [iTOL](https://itol.embl.de/)
2. Upload `results/tree.treefile`
3. Customize colors, labels, scale

**Locally:**
```bash
# Install FigTree (macOS)
brew install figtree
figtree results/tree.treefile
```

## Output Files

### Balanced Pipeline Output

| File | Description |
|------|-------------|
| `sequences_balanced.fasta` | Equal sequences per species (~5 each) |
| `sequences_clustered.fasta` | After CD-HIT redundancy removal |
| `sequences_aligned.fasta` | MAFFT multiple sequence alignment |
| `tree.treefile` | Phylogenetic tree (Newick format) |
| `pipeline.log` | Full execution log |

### Tree-Guided Pipeline Output

| File | Description |
|------|-------------|
| `sequences_nonredundant.fasta` | After CD-HIT clustering |
| `sequences_aligned.fasta` | MAFFT multiple sequence alignment |
| `diverse_dataset.fasta` | Final selected sequences |
| `tree.treefile` | Phylogenetic tree |
| `pipeline.log` | Full execution log |

## Configuration

### Download Custom Organisms

Edit `config/config.json`:

```json
{
  "species": {
    "Species1": "organism_name:Exact species name",
    "Species2": "organism_name:Another species"
  },
  "protein_filter": "protein_name:envelope glycoprotein"
}
```

Then download:
```bash
docker compose run ebolavirus-analysis \
  python scripts/download_gp_uniprot.py \
  -o sequences \
  -c config/config.json
```

### Pipeline Command-Line Options

#### Balanced Diversity Pipeline

```bash
python scripts/balanced_diversity_pipeline.py \
  -i sequences.fasta \
  -o results \
  --target-per-species 5 \
  --cd-hit-threshold 0.80 \
  --seq-type protein
```

Options:
- `-i, --input` - Input FASTA file (required)
- `-o, --output` - Output directory (default: results)
- `--target-per-species` - Sequences per species (default: auto)
- `--cd-hit-threshold` - Identity threshold 0.60-0.95 (default: 0.80)
- `--seq-type` - `protein` or `nucleotide` (default: protein)

#### Tree-Guided Pipeline

```bash
python scripts/generic_pipeline.py \
  -i sequences.fasta \
  -o results \
  -n 40 \
  --cd-hit-threshold 0.95 \
  --use-fasttree \
  --seq-type protein
```

Options:
- `-i, --input` - Input FASTA file (required)
- `-o, --output` - Output directory (default: results)
- `-n, --num-sequences` - Target sequences to select (default: 40)
- `--cd-hit-threshold` - Identity threshold (default: 0.95)
- `--use-fasttree` - Use FastTree instead of IQ-TREE (faster)
- `--seq-type` - `protein` or `nucleotide` (default: protein)

## Advanced Usage

### Run Individual Steps

```bash
# Just clustering
docker compose run ebolavirus-analysis bash -c \
  "cd-hit -i sequences/ebolavirus_gp_combined.fasta \
   -o results/nr.fasta -c 0.95 -n 5"

# Just alignment
docker compose run ebolavirus-analysis bash -c \
  "mafft --auto results/nr.fasta > results/aligned.fasta"

# Just tree building (IQ-TREE)
docker compose run ebolavirus-analysis bash -c \
  "iqtree -s results/aligned.fasta -st AA -m LG+G -bb 1000 -nt AUTO"
```

### Use FastTree for Speed

For large datasets, FastTree is 10-100x faster than IQ-TREE:

```bash
docker compose run ebolavirus-analysis \
  python scripts/generic_pipeline.py \
  -i sequences.fasta \
  -n 40 \
  --use-fasttree
```

### Interactive Shell

For manual analysis:

```bash
docker compose run ebolavirus-analysis bash
```

Inside container:
```bash
cd /data
mafft --version
iqtree --version
cd-hit --version
```

### Nucleotide Sequences

Both pipelines support DNA/RNA:

```bash
# Balanced selection for nucleotides
python scripts/balanced_diversity_pipeline.py \
  -i dna_sequences.fasta \
  --seq-type nucleotide

# Tree-guided for nucleotides
python scripts/generic_pipeline.py \
  -i dna_sequences.fasta \
  --seq-type nucleotide
```

## When to Use Which Pipeline

| Scenario | Use |
|----------|-----|
| Imbalanced species (e.g., 400 Zaire, 3 Bundibugyo) | **Balanced** |
| Large dataset, want maximum diversity | **Tree-guided** |
| Already balanced species representation | Either |
| Small dataset (<50 seqs) | **Balanced** |
| Massive dataset (>1000 seqs) | **Tree-guided** with FastTree |

## Understanding CD-HIT Thresholds

- **0.95 (default for tree-guided)**: Very aggressive, removes more sequences
- **0.80 (default for balanced)**: Moderate, balances diversity vs redundancy removal
- **0.70**: Less aggressive, keeps more sequences but may retain near-duplicates
- **0.60**: Minimal clustering, keep almost everything

**Lower threshold = more sequences kept = slower analysis but better diversity**

## Troubleshooting

### Out of Memory

Increase Docker memory:
- Docker Desktop: Settings → Resources → Memory → 8GB+
- Or edit `docker-compose.yml`:
  ```yaml
  services:
    ebolavirus-analysis:
      mem_limit: 8g
  ```

### Pipeline Takes Too Long

Use FastTree instead of IQ-TREE:
```bash
python scripts/generic_pipeline.py -i sequences.fasta --use-fasttree
```

### Missing Sequences in Results

Check the log file:
```bash
docker compose run ebolavirus-analysis bash -c \
  "tail -100 results/pipeline.log"
```

Common issues:
- Sequences too similar (try lower `--cd-hit-threshold`)
- Corrupted FASTA file (verify with `docker compose run ebolavirus-analysis bash -c "head -20 sequences/your_file.fasta"`)
- Alignment failed (check sequence lengths vary too much)

### CD-HIT Not Removing Redundancy

Lower the threshold - 0.95 may be too high for your sequences:
```bash
python scripts/balanced_diversity_pipeline.py \
  -i sequences.fasta \
  --cd-hit-threshold 0.80
```

### "docker compose: command not found"

Update Docker or use older syntax:
```bash
docker-compose run ebolavirus-analysis ...  # Old
docker compose run ebolavirus-analysis ...   # New (preferred)
```

## Clean Up

```bash
# Remove old containers
docker compose down --remove-orphans

# Delete results and start over
rm -rf results/*

# Clean everything
docker compose down -v
rm -rf results sequences
```

## Next Steps with Your Data

Once you have your selected sequences:

1. **3D Structure Prediction** - [ColabFold](https://github.com/sokrypton/ColabFold)
2. **Sequence Conservation** - [ConSurf](https://consurf.tau.ac.il/)
3. **Epitope Mapping** - IEDB, BepiPred
4. **Molecular Dynamics** - GROMACS, NAMD
5. **Domain Analysis** - InterProScan, Pfam
6. **Comparative Analysis** - Multiple sequence comparison tools

## For Ebolavirus Specifically

The original ebolavirus-specific scripts still work:

```bash
# Download ebolavirus sequences
docker compose run ebolavirus-analysis \
  python scripts/download_gp_uniprot.py -o sequences

# Analyze species coverage
docker compose run ebolavirus-analysis \
  python scripts/analyze_subtypes.py sequences/ebolavirus_gp_combined.fasta

# Original pipeline (less flexible but still works)
docker compose run ebolavirus-analysis \
  python scripts/pipeline.py \
  -i sequences/ebolavirus_gp_combined.fasta \
  -o results \
  -n 40
```

## Citations

If you use this pipeline, please cite:

- **MAFFT**: Katoh K, et al. (2005) MAFFT: a novel method for rapid multiple sequence alignment
- **IQ-TREE**: Nguyen LT, Schmidt HA, et al. (2015) IQ-TREE: A fast and effective stochastic algorithm
- **FastTree**: Price MN, et al. (2010) FastTree 2: approximately maximum-likelihood trees for large alignments
- **CD-HIT**: Li W, Godzik A. (2006) Cd-hit: a fast program for clustering and comparing large sets
- **Biopython**: Cock PJ, et al. (2009) Biopython: freely available Python tools for computational molecular biology

## Support

For issues:

1. Check `results/pipeline.log` for error details
2. Review troubleshooting section above
3. Verify Docker: `docker compose version`
4. Check disk/memory: `df -h` and `docker stats`
5. Test internet: `docker compose run ebolavirus-analysis bash -c "curl -I https://rest.uniprot.org/"`

## License

This pipeline is provided for research purposes. See individual tool licenses for terms.

---

**Happy analyzing!** 🧬
