# Ebolavirus Docker Container - Complete Tutorial

This guide walks you through using the Docker container from start to finish.

## System Requirements

- **Docker Desktop** (or Docker + Docker Compose on Linux)
- **4GB RAM** minimum (8GB recommended)
- **20GB disk space** (for container + sequences + results)
- **Internet connection** (for initial build and downloading sequences)

## Installation

### Step 1: Install Docker

**macOS/Windows:**
- Download and install [Docker Desktop](https://www.docker.com/products/docker-desktop)

**Linux (Ubuntu/Debian):**
```bash
sudo apt-get update
sudo apt-get install docker.io docker-compose
sudo usermod -aG docker $USER
# Logout and login again for group changes to take effect
```

### Step 2: Clone or setup the project

```bash
# Create project directory
mkdir ebolavirus-analysis
cd ebolavirus-analysis

# Copy these files (from earlier):
# - Dockerfile
# - docker-compose.yml
# - pipeline.py
# - utils.py
# - entrypoint.sh
# - setup.sh
# - README.md

# Make setup script executable
chmod +x setup.sh

# Run setup
./setup.sh
```

This will:
- Create necessary directories (sequences/, results/, scripts/)
- Build the Docker image (takes ~10-15 minutes the first time)
- Prepare the environment

## Complete Example Workflow

### Example 1: Analyze Downloaded Sequences

Suppose you've downloaded ebolavirus sequences from NCBI as `ebolavirus.fasta`.

**Step 1: Place sequences in container directory**

```bash
cp ebolavirus.fasta sequences/
```

**Step 2: Run the analysis pipeline**

```bash
docker-compose run ebolavirus-analysis \
  python scripts/pipeline.py \
  -i sequences/ebolavirus.fasta \
  -o results \
  -n 40
```

This will:
1. Remove redundant sequences (CD-HIT at 95% identity)
2. Align them (MAFFT)
3. Build a phylogenetic tree (IQ-TREE)
4. Select 40 diverse sequences guided by the tree
5. Save results to `results/`

**Step 3: Check the output**

```bash
# View what was generated
ls -lh results/

# See how many sequences were selected
grep "^>" results/ebolavirus_diverse_dataset.fasta | wc -l

# View tree structure
cat results/tree.treefile

# Check pipeline log
cat results/pipeline.log | tail -20
```

### Example 2: Interactive Analysis

If you want more control, use interactive bash:

```bash
docker-compose run ebolavirus-analysis bash
```

Now you're inside the container. You have access to all tools:

```bash
# Inside container:

# Check MAFFT version
mafft --version

# Manually run CD-HIT on your sequences
cd-hit -i /data/sequences/ebolavirus.fasta \
        -o /data/sequences/ebolavirus_nr.fasta \
        -c 0.95 -n 8

# Manually align
mafft --auto /data/sequences/ebolavirus_nr.fasta > /data/sequences/aligned.fasta

# Build tree with IQ-TREE
iqtree -s /data/sequences/aligned.fasta -m GTR+G -bb 1000

# Exit container
exit
```

### Example 3: Download Sequences Inside Container

Instead of downloading locally, do it inside the container:

```bash
docker-compose run ebolavirus-analysis bash

# Inside container:
# Install NCBI datasets tool
cd /data

# Download ebolavirus sequences
datasets download virus genome taxon ebolavirus --assembly-type all

# Extract FASTA
unzip ncbi_dataset.zip
cat ncbi_dataset/data/genomic.fasta > sequences/ebolavirus.fasta

# Run pipeline
python scripts/pipeline.py -i sequences/ebolavirus.fasta -o results -n 40

exit
```

### Example 4: Analyze Results with Utilities

```bash
# Analyze the alignment
docker-compose run ebolavirus-analysis \
  python scripts/utils.py --alignment results/sequences_aligned.fasta

# Analyze the tree
docker-compose run ebolavirus-analysis \
  python scripts/utils.py --tree results/tree.treefile

# Get sequence statistics
docker-compose run ebolavirus-analysis \
  python scripts/utils.py --stats results/ebolavirus_diverse_dataset.fasta

# Generate full report
docker-compose run ebolavirus-analysis \
  python scripts/utils.py --report results/
```

## Common Workflows

### Workflow A: Quick Analysis

```bash
# 1. Download sequences locally
# 2. Copy to sequences/
cp my_sequences.fasta sequences/

# 3. Run one-liner
docker-compose run ebolavirus-analysis \
  python scripts/pipeline.py \
  -i sequences/my_sequences.fasta \
  -o results -n 50

# 4. Results in results/ directory on your machine
```

### Workflow B: Large Dataset (1000+ sequences)

```bash
# Use FastTree for speed
docker-compose run ebolavirus-analysis \
  python scripts/pipeline.py \
  -i sequences/large_dataset.fasta \
  --use-fasttree \
  -o results \
  -n 100
```

### Workflow C: Manual Step-by-Step Control

```bash
docker-compose run ebolavirus-analysis bash

# Inside container:
cd /data

# Step 1: Redundancy removal
cd-hit -i sequences/raw.fasta \
        -o sequences/nonredundant.fasta \
        -c 0.98 -n 8

# Step 2: Alignment
mafft --auto sequences/nonredundant.fasta > results/aligned.fasta

# Step 3: Inspect alignment
seqmagick info results/aligned.fasta

# Step 4: Tree building
iqtree -s results/aligned.fasta -m GTR+G -bb 1000 -nt 4

# Step 5: Analyze with Python
python << 'EOF'
from Bio import SeqIO, Phylo
tree = Phylo.read("results/aligned.treefile", "newick")
Phylo.draw_ascii(tree)
EOF

exit
```

## Understanding the Output Files

After running the pipeline, you'll have:

```
results/
├── sequences_nonredundant.fasta
│   └── All sequences after removing 95%+ identical ones
│
├── sequences_aligned.fasta
│   └── Multiple sequence alignment (can view in MEGA, AliView, etc.)
│
├── tree.treefile
│   └── Phylogenetic tree (Newick format)
│       Can open with: FigTree, iTOL, or R
│
├── ebolavirus_diverse_dataset.fasta
│   └── YOUR FINAL RESULT: 40 (or -n value) diverse representatives
│
├── tree.iqtree
│   └── IQ-TREE detailed output
│
└── pipeline.log
    └── Execution log with all steps and statistics
```

## Visualizing Results

### View the tree online

1. Go to [iTOL](https://itol.embl.de/)
2. Upload `results/tree.treefile`
3. Customize visualization (colors, labels, etc.)

### View the tree locally

Download and install [FigTree](http://tree.bio.ed.ac.uk/software/figtree/):

```bash
# macOS
brew install figtree

# Then open the tree file
figtree results/tree.treefile
```

### View the alignment

Download [AliView](http://www.ormbunkar.se/aliview/):

```bash
# Open alignment
aliview results/sequences_aligned.fasta
```

## Troubleshooting

### Problem: "docker-compose: command not found"

**Solution (macOS/Windows):** Docker Desktop includes docker-compose  
**Solution (Linux):** Install it

```bash
sudo apt-get install docker-compose
# Or
pip install docker-compose
```

### Problem: "Permission denied while trying to connect to Docker daemon"

**Solution (Linux):**

```bash
sudo usermod -aG docker $USER
# Then logout and login again
```

### Problem: "No such file or directory: sequences/ebolavirus.fasta"

**Solution:** Make sure you've placed the FASTA file in the `sequences/` directory:

```bash
ls -la sequences/
# Should show your FASTA file
```

### Problem: "out of memory" errors

**Solution:** Increase Docker's memory allocation

**Docker Desktop (GUI):**
- Settings → Resources → Memory → Increase to 8GB+

**Or edit docker-compose.yml:**

```yaml
services:
  ebolavirus-analysis:
    mem_limit: 8g  # Add this line
```

### Problem: Tree building takes forever

**Solution:** Use FastTree instead of IQ-TREE:

```bash
docker-compose run ebolavirus-analysis \
  python scripts/pipeline.py \
  --use-fasttree \
  -i sequences/ebolavirus.fasta \
  -o results
```

## Customization Examples

### Select a different number of sequences

```bash
# Get 100 diverse sequences instead of 40
docker-compose run ebolavirus-analysis \
  python scripts/pipeline.py \
  -i sequences/ebolavirus.fasta \
  -o results \
  -n 100
```

### Use a different clustering threshold

Edit `scripts/pipeline.py`, find this line:

```python
self.remove_redundancy(input_fasta, str(nonredundant), identity_threshold=0.95)
```

Change `0.95` to:
- `0.99` - More aggressive, removes even minor variants
- `0.90` - Less aggressive, keeps more diversity

### Add custom Python analysis

Create `scripts/my_analysis.py`:

```python
#!/usr/bin/env python3
from Bio import SeqIO, Phylo

# Analyze the final dataset
for record in SeqIO.parse("/data/results/ebolavirus_diverse_dataset.fasta", "fasta"):
    print(f"{record.id}: {len(record.seq)} bp")

# Analyze the tree
tree = Phylo.read("/data/results/tree.treefile", "newick")
terminals = [c for c in tree.find_clades(terminal=True)]
print(f"Tree has {len(terminals)} sequences")
```

Run it:

```bash
docker-compose run ebolavirus-analysis python scripts/my_analysis.py
```

## Data Management

### Backup your results

```bash
tar -czf ebolavirus_results_$(date +%Y%m%d).tar.gz results/
```

### Clean up old results

```bash
# Keep sequences/, but remove old results
rm -rf results/*

# But first backup:
cp -r results results_backup_$(date +%Y%m%d)
```

### Share results with collaborators

```bash
# Share the alignment
scp results/sequences_aligned.fasta colleague@server:/shared/

# Share the tree
scp results/tree.treefile colleague@server:/shared/

# Share the final dataset
scp results/ebolavirus_diverse_dataset.fasta colleague@server:/shared/
```

## Performance Notes

### Typical runtimes on modern machines

| Step | Time | Tool |
|------|------|------|
| CD-HIT (redundancy) | 1-2 min | CD-HIT |
| Alignment (100 seqs) | 2-5 min | MAFFT |
| Tree (IQ-TREE) | 10-30 min | IQ-TREE |
| Tree (FastTree) | 1-2 min | FastTree |
| Phylo-selection | <1 min | Python |
| **Total (default)** | **20-40 min** | - |
| **Total (FastTree)** | **5-10 min** | - |

### Memory usage

- Typical: 2-4 GB RAM
- Large datasets (1000+ seqs): 6-8 GB RAM

## Next Steps

After getting your diverse dataset:

1. **Sequence alignment analysis** - Use MEGA or SeaView for detailed conservation analysis
2. **Molecular clock analysis** - Estimate divergence times with BEAST or r8s
3. **Recombination detection** - Use RDP4 or Simplot
4. **Protein structure prediction** - Use ColabFold with your aligned sequences
5. **Comparative genomics** - Analyze synteny, GC content, codon usage

## Citing the Pipeline

If you use this pipeline in publications, please cite:

- **MAFFT**: Katoh K, Misawa K, et al. (2002) MAFFT: a novel method for rapid multiple sequence alignment
- **IQ-TREE**: Nguyen LT, Schmidt HA, et al. (2015) IQ-TREE: fast and effective stochastic algorithm
- **CD-HIT**: Li W, Godzik A. (2006) Cd-hit: a fast program for clustering
- **Biopython**: Cock PJ, et al. (2009) Biopython: freely available Python tools

## Support & Feedback

- Check `results/pipeline.log` for detailed execution information
- Review the README.md for comprehensive documentation
- All tools (MAFFT, IQ-TREE, etc.) have their own documentation

---

Happy analyzing! 🔬
