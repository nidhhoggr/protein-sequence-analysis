#!/usr/bin/env python3
"""
Balanced Species Diversity Selection Pipeline

Creates a phylogenetically-informed, species-balanced dataset for any organism.
Ensures equal representation across all species in your dataset.

Solves the problem: When you have many sequences from one species and few from others,
this pipeline selects equal numbers from each species for unbiased phylogenetic analysis.

Usage:
  python balanced_diversity_pipeline.py -i sequences.fasta -o results
  
  Optional:
    --cd-hit-threshold 0.80   (default: 0.80, range: 0.60-0.95)
    --target-per-species 5    (default: auto-calculate from -n)
    --seq-type protein        (default: protein, or 'nucleotide')
"""

import os
import sys
import subprocess
import argparse
import re
from pathlib import Path
from Bio import SeqIO, Phylo
from collections import defaultdict

class BalancedDiversityAnalysis:
    def __init__(self, output_dir="results", num_sequences=None, cd_hit_threshold=0.80, 
                 target_per_species=None, seq_type="protein"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.num_sequences = num_sequences
        self.cd_hit_threshold = cd_hit_threshold
        self.target_per_species = target_per_species
        self.seq_type = seq_type  # 'protein' or 'nucleotide'
        self.log_file = self.output_dir / "pipeline.log"
        
    def log(self, message):
        """Log to file and stdout"""
        print(message)
        with open(self.log_file, "a") as f:
            f.write(message + "\n")
    
    def run_command(self, cmd, description):
        """Execute shell command"""
        self.log(f"\n{'='*60}")
        self.log(f"Running: {description}")
        self.log(f"Command: {cmd}")
        self.log(f"{'='*60}")
        
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode != 0:
                self.log(f"ERROR: {result.stderr}")
                return False
            # Log last 500 chars to avoid spam
            output = result.stdout[-500:] if len(result.stdout) > 500 else result.stdout
            self.log(output)
            return True
        except Exception as e:
            self.log(f"ERROR: {str(e)}")
            return False
    
    def extract_species(self, header):
        """
        Extract species from FASTA header.
        Looks for 'OS=' pattern (UniProt format).
        Falls back to looking for pattern like "Species_name" or just uses full header.
        """
        # Try UniProt format: OS=Species name
        match = re.search(r'OS=([^O]+?)(?:\s+OX=|$)', header)
        if match:
            sp = match.group(1).strip()
            # Remove strain info in parentheses
            if '(strain' in sp or '(isolate' in sp:
                sp = re.sub(r'\s*\(.*\)\s*', '', sp).strip()
            return sp
        
        # Fallback: try to extract something meaningful
        # Look for pattern like ">species_name" or first part before pipe
        parts = header.split('|')
        if len(parts) > 1:
            return parts[1].strip()
        
        return "Unknown"
    
    def select_balanced_sequences(self, input_fasta, output_fasta):
        """Select balanced representatives from each species"""
        self.log(f"\nAnalyzing sequence diversity...")
        
        # Group by species
        by_species = defaultdict(list)
        total_seqs = 0
        
        for record in SeqIO.parse(input_fasta, 'fasta'):
            sp = self.extract_species(record.description)
            by_species[sp].append(record)
            total_seqs += 1
        
        num_species = len(by_species)
        self.log(f"Found {num_species} species in {total_seqs} sequences")
        self.log(f"\nSpecies breakdown:")
        
        for sp in sorted(by_species.keys()):
            count = len(by_species[sp])
            pct = 100 * count / total_seqs
            self.log(f"  {sp:40s}: {count:4d} ({pct:5.1f}%)")
        
        # Determine target per species
        if self.target_per_species is None:
            if self.num_sequences:
                self.target_per_species = max(1, self.num_sequences // num_species)
            else:
                self.target_per_species = 5  # Default
        
        self.log(f"\nSelecting {self.target_per_species} sequences per species:")
        
        # Select evenly-spaced from each species
        balanced_seqs = []
        for sp in sorted(by_species.keys()):
            seqs = by_species[sp]
            # Select evenly-spaced to maximize diversity within species
            step = max(1, len(seqs) // self.target_per_species)
            selected = seqs[::step][:self.target_per_species]
            balanced_seqs.extend(selected)
            self.log(f"  {sp:40s}: {len(selected):2d} of {len(seqs):4d}")
        
        self.log(f"\nTotal selected: {len(balanced_seqs)} sequences")
        SeqIO.write(balanced_seqs, output_fasta, 'fasta')
        
        return len(balanced_seqs) > 0
    
    def cluster_sequences(self, input_fasta, output_fasta):
        """CD-HIT clustering to remove near-duplicates"""
        self.log(f"\nRemoving near-duplicates with CD-HIT at {self.cd_hit_threshold*100}% identity...")
        
        # CD-HIT word length: 5 for proteins, 10 for nucleotides
        word_length = 5 if self.seq_type == "protein" else 10
        
        cmd = f"cd-hit -i {input_fasta} -o {output_fasta} -c {self.cd_hit_threshold} -n {word_length} -d 0 -T 0 2>&1"
        
        if not self.run_command(cmd, "CD-HIT clustering"):
            return False
        
        # Count output
        cluster_count = sum(1 for _ in SeqIO.parse(output_fasta, 'fasta'))
        self.log(f"Result: {cluster_count} representative sequences after clustering")
        
        return True
    
    def align_sequences(self, input_fasta, output_fasta):
        """MAFFT multiple sequence alignment"""
        self.log(f"\nAligning sequences with MAFFT...")
        
        cmd = f"mafft --auto {input_fasta} > {output_fasta} 2>&1"
        return self.run_command(cmd, "MAFFT alignment")
    
    def build_tree(self, aligned_fasta, output_prefix):
        """IQ-TREE phylogenetic inference"""
        self.log(f"\nBuilding phylogenetic tree with IQ-TREE...")
        
        # Use protein model for amino acids, nucleotide for DNA/RNA
        seq_type_flag = "AA" if self.seq_type == "protein" else "DNA"
        
        cmd = f"iqtree -s {aligned_fasta} -st {seq_type_flag} -m LG+G -bb 1000 -nt AUTO -pre {output_prefix} 2>&1"
        return self.run_command(cmd, "IQ-TREE tree construction")
    
    def run_pipeline(self, input_fasta):
        """Execute full balanced diversity pipeline"""
        self.log(f"Starting Balanced Diversity Analysis Pipeline")
        self.log(f"{'='*60}")
        self.log(f"Input: {input_fasta}")
        self.log(f"Output: {self.output_dir}")
        self.log(f"Sequence type: {self.seq_type}")
        self.log(f"Target per species: {self.target_per_species or 'auto'}")
        self.log(f"CD-HIT threshold: {self.cd_hit_threshold*100}%")
        
        # Verify input exists
        if not Path(input_fasta).exists():
            self.log(f"ERROR: Input file not found: {input_fasta}")
            return False
        
        # Step 1: Balance by species
        balanced = self.output_dir / "sequences_balanced.fasta"
        if not self.select_balanced_sequences(input_fasta, str(balanced)):
            return False
        
        # Step 2: Cluster to remove duplicates
        clustered = self.output_dir / "sequences_clustered.fasta"
        if not self.cluster_sequences(str(balanced), str(clustered)):
            return False
        
        # Step 3: Align
        aligned = self.output_dir / "sequences_aligned.fasta"
        if not self.align_sequences(str(clustered), str(aligned)):
            return False
        
        # Step 4: Build tree
        tree_prefix = str(self.output_dir / "tree")
        if not self.build_tree(str(aligned), tree_prefix):
            return False
        
        self.log(f"\n{'='*60}")
        self.log("Pipeline Complete!")
        self.log(f"{'='*60}")
        self.log(f"\nOutput files:")
        self.log(f"  Balanced dataset: {balanced}")
        self.log(f"  Clustered: {clustered}")
        self.log(f"  Alignment: {aligned}")
        self.log(f"  Phylogenetic tree: {tree_prefix}.treefile")
        self.log(f"  Full log: {self.log_file}")
        
        return True


def main():
    parser = argparse.ArgumentParser(
        description="Balanced Species Diversity Selection Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic: auto-balance all species
  python balanced_diversity_pipeline.py -i sequences.fasta
  
  # Custom: 7 sequences per species
  python balanced_diversity_pipeline.py -i sequences.fasta --target-per-species 7
  
  # Lower CD-HIT threshold for more diversity
  python balanced_diversity_pipeline.py -i sequences.fasta --cd-hit-threshold 0.70
  
  # For nucleotide sequences
  python balanced_diversity_pipeline.py -i dna_sequences.fasta --seq-type nucleotide
  
  # Custom output directory
  python balanced_diversity_pipeline.py -i sequences.fasta -o my_results
        """
    )
    
    parser.add_argument("-i", "--input", required=True,
                        help="Input FASTA file")
    parser.add_argument("-o", "--output", default="results",
                        help="Output directory (default: results)")
    parser.add_argument("-n", "--num-sequences", type=int, default=None,
                        help="Target total sequences (optional)")
    parser.add_argument("--cd-hit-threshold", type=float, default=0.80,
                        help="CD-HIT identity threshold (default: 0.80, range: 0.60-0.95)")
    parser.add_argument("--target-per-species", type=int, default=None,
                        help="Sequences per species (default: auto-calculated)")
    parser.add_argument("--seq-type", choices=["protein", "nucleotide"], default="protein",
                        help="Sequence type (default: protein)")
    
    args = parser.parse_args()
    
    pipeline = BalancedDiversityAnalysis(
        output_dir=args.output,
        num_sequences=args.num_sequences,
        cd_hit_threshold=args.cd_hit_threshold,
        target_per_species=args.target_per_species,
        seq_type=args.seq_type
    )
    
    success = pipeline.run_pipeline(args.input)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
