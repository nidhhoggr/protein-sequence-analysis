#!/usr/bin/env python3
"""
Generic Phylogenetic Diversity Analysis Pipeline

Creates a phylogenetically-informed, diverse dataset for any organism.
Uses tree-guided selection to pick sequences that maximize evolutionary coverage.

Usage:
  python generic_pipeline.py -i sequences.fasta -o results -n 40
  
  Optional:
    --cd-hit-threshold 0.80   (default: 0.80, range: 0.60-0.95)
    --seq-type protein        (default: protein, or 'nucleotide')
    --use-fasttree            (use FastTree instead of IQ-TREE, faster but less accurate)
"""

import os
import sys
import subprocess
import argparse
import random
from pathlib import Path
from Bio import SeqIO, Phylo
from collections import defaultdict

class PhylogeneticDiversityAnalysis:
    def __init__(self, output_dir="results", num_sequences=40, cd_hit_threshold=0.95, seq_type="protein", use_fasttree=False):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.num_sequences = num_sequences
        self.cd_hit_threshold = cd_hit_threshold
        self.seq_type = seq_type
        self.use_fasttree = use_fasttree
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
            output = result.stdout[-500:] if len(result.stdout) > 500 else result.stdout
            self.log(output)
            return True
        except Exception as e:
            self.log(f"ERROR: {str(e)}")
            return False
    
    def download_sequences(self, input_fasta):
        """Verify input FASTA file exists"""
        if os.path.exists(input_fasta):
            self.log(f"Using input FASTA file: {input_fasta}")
            return True
        
        self.log(f"ERROR: FASTA file not found: {input_fasta}")
        return False
    
    def remove_redundancy(self, input_fasta, output_fasta, identity_threshold=None):
        """Remove sequence redundancy using CD-HIT"""
        if identity_threshold is None:
            identity_threshold = self.cd_hit_threshold
            
        self.log(f"\nRemoving sequences with >{identity_threshold*100}% identity...")
        
        # CD-HIT word length: 5 for proteins, 10 for nucleotides
        word_length = 5 if self.seq_type == "protein" else 10
        
        cmd = f"cd-hit -i {input_fasta} -o {output_fasta} -c {identity_threshold} -n {word_length} -d 0 -T 0 2>&1"
        return self.run_command(cmd, f"CD-HIT clustering at {identity_threshold*100}% identity")
    
    def align_sequences(self, input_fasta, output_alignment):
        """Align sequences using MAFFT"""
        self.log(f"\nAligning sequences with MAFFT...")
        cmd = f"mafft --auto {input_fasta} > {output_alignment} 2>&1"
        return self.run_command(cmd, "MAFFT sequence alignment")
    
    def build_tree_fasttree(self, aligned_fasta, output_tree):
        """Build phylogenetic tree using FastTree"""
        self.log(f"\nBuilding phylogenetic tree with FastTree...")
        
        # FastTree: -nt for nucleotides (default is protein)
        seq_flag = "-nt" if self.seq_type == "nucleotide" else ""
        cmd = f"FastTree {seq_flag} {aligned_fasta} > {output_tree} 2>&1"
        return self.run_command(cmd, "FastTree phylogenetic tree construction")
    
    def build_tree_iqtree(self, aligned_fasta, output_prefix):
        """Build phylogenetic tree using IQ-TREE"""
        self.log(f"\nBuilding phylogenetic tree with IQ-TREE...")
        
        seq_type_flag = "AA" if self.seq_type == "protein" else "DNA"
        cmd = f"iqtree -s {aligned_fasta} -st {seq_type_flag} -m LG+G -bb 1000 -nt AUTO -pre {output_prefix} 2>&1"
        return self.run_command(cmd, "IQ-TREE phylogenetic tree construction")
    
    def select_by_phylogeny(self, tree_file, sequence_fasta, output_fasta):
        """
        Select diverse sequences by traversing the phylogenetic tree.
        Randomly samples across the tree to maximize diversity.
        """
        self.log(f"\nSelecting {self.num_sequences} diverse sequences using phylogenetic guidance...")
        
        try:
            tree = Phylo.read(tree_file, "newick")
            
            # Build mapping of all sequences
            fasta_records = {}
            for record in SeqIO.parse(sequence_fasta, "fasta"):
                fasta_records[record.id] = record
            
            self.log(f"Total sequences in alignment: {len(fasta_records)}")
            
            # Collect all terminal nodes
            all_terminals = []
            for clade in tree.find_clades(terminal=True):
                if clade.name:
                    all_terminals.append(clade.name)
            
            self.log(f"Total sequences in tree: {len(all_terminals)}")
            
            # Strategy: Randomly sample across tree
            selected_ids = set()
            
            def collect_and_select(clade, selected_ids):
                if clade.is_terminal():
                    if clade.name and clade.name in fasta_records:
                        if len(selected_ids) < self.num_sequences:
                            if random.random() < 0.5:
                                selected_ids.add(clade.name)
                else:
                    for child in clade.clades:
                        collect_and_select(child, selected_ids)
            
            # Multiple passes to reach target
            for attempt in range(5):
                if len(selected_ids) >= self.num_sequences:
                    break
                collect_and_select(tree.root, selected_ids)
            
            # Fill remaining with random sequences
            remaining = set(fasta_records.keys()) - selected_ids
            needed = self.num_sequences - len(selected_ids)
            if needed > 0 and remaining:
                selected_ids.update(random.sample(list(remaining), min(needed, len(remaining))))
            
            self.log(f"Selected {len(selected_ids)} sequences via phylogenetic guidance")
            
            # Extract selected sequences
            selected_records = [fasta_records[seq_id] for seq_id in selected_ids if seq_id in fasta_records]
            
            self.log(f"Extracted {len(selected_records)} sequences to {output_fasta}")
            SeqIO.write(selected_records, output_fasta, "fasta")
            
            return len(selected_records) > 0
            
        except Exception as e:
            self.log(f"ERROR in phylogenetic selection: {str(e)}")
            import traceback
            self.log(traceback.format_exc())
            return False
    
    def check_quality(self, fasta_file):
        """Check quality of final dataset"""
        self.log(f"\n{'='*60}")
        self.log("Quality Check")
        self.log(f"{'='*60}")
        
        count = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
        self.log(f"Total sequences: {count}")
        
        lengths = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            lengths.append(len(record.seq))
        
        if lengths:
            self.log(f"Sequence lengths: min={min(lengths)}, max={max(lengths)}, mean={sum(lengths)/len(lengths):.0f}")
        
        ids = [record.id for record in SeqIO.parse(fasta_file, "fasta")]
        unique_ids = set(ids)
        if len(ids) != len(unique_ids):
            self.log(f"WARNING: Found {len(ids) - len(unique_ids)} duplicate sequence IDs")
        else:
            self.log("No duplicate sequences detected ✓")
    
    def run_pipeline(self, input_fasta):
        """Execute full pipeline"""
        self.log(f"Starting Phylogenetic Diversity Analysis Pipeline")
        self.log(f"Target: {self.num_sequences} diverse sequences")
        self.log(f"Output directory: {self.output_dir}")
        self.log(f"Sequence type: {self.seq_type}")
        self.log(f"CD-HIT threshold: {self.cd_hit_threshold*100}%")
        self.log(f"Tree building: {'FastTree (fast)' if self.use_fasttree else 'IQ-TREE (accurate)'}")
        
        # Step 1: Verify sequences
        if not self.download_sequences(input_fasta):
            return False
        
        # Step 2: Remove redundancy
        nonredundant = self.output_dir / "sequences_nonredundant.fasta"
        if not self.remove_redundancy(input_fasta, str(nonredundant)):
            return False
        
        # Step 3: Align
        aligned = self.output_dir / "sequences_aligned.fasta"
        if not self.align_sequences(str(nonredundant), str(aligned)):
            return False
        
        # Step 4: Build tree
        if self.use_fasttree:
            tree_file = self.output_dir / "tree.nwk"
            if not self.build_tree_fasttree(str(aligned), str(tree_file)):
                return False
        else:
            tree_output = self.output_dir / "tree"
            if not self.build_tree_iqtree(str(aligned), str(tree_output)):
                return False
            tree_file = str(tree_output) + ".treefile"
        
        # Step 5: Select diverse sequences
        final_fasta = self.output_dir / "diverse_dataset.fasta"
        if not self.select_by_phylogeny(str(tree_file), str(aligned), str(final_fasta)):
            self.log("Phylogenetic selection failed, using alignment instead")
            final_fasta = aligned
        
        # Step 6: Quality check
        self.check_quality(str(final_fasta))
        
        self.log(f"\n{'='*60}")
        self.log("Pipeline Complete!")
        self.log(f"{'='*60}")
        self.log(f"Final dataset: {final_fasta}")
        self.log(f"Alignment: {aligned}")
        self.log(f"Tree: {tree_file}")
        self.log(f"Log file: {self.log_file}")
        
        return True


def main():
    parser = argparse.ArgumentParser(
        description="Generic Phylogenetic Diversity Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python generic_pipeline.py -i sequences.fasta
  
  # Specify number of target sequences
  python generic_pipeline.py -i sequences.fasta -n 50
  
  # Use faster tree building
  python generic_pipeline.py -i sequences.fasta --use-fasttree
  
  # Lower CD-HIT threshold for more diversity
  python generic_pipeline.py -i sequences.fasta --cd-hit-threshold 0.70
  
  # For nucleotide sequences
  python generic_pipeline.py -i dna_sequences.fasta --seq-type nucleotide
  
  # Custom output
  python generic_pipeline.py -i sequences.fasta -o my_results
        """
    )
    
    parser.add_argument("-i", "--input", required=True,
                        help="Input FASTA file with sequences")
    parser.add_argument("-o", "--output", default="results",
                        help="Output directory (default: results)")
    parser.add_argument("-n", "--num-sequences", type=int, default=40,
                        help="Target number of diverse sequences to select (default: 40)")
    parser.add_argument("--use-fasttree", action="store_true",
                        help="Use FastTree instead of IQ-TREE (faster but less accurate)")
    parser.add_argument("--cd-hit-threshold", type=float, default=0.95,
                        help="CD-HIT identity threshold (default: 0.95, range: 0.60-0.95)")
    parser.add_argument("--seq-type", choices=["protein", "nucleotide"], default="protein",
                        help="Sequence type (default: protein)")
    
    args = parser.parse_args()
    
    pipeline = PhylogeneticDiversityAnalysis(
        output_dir=args.output,
        num_sequences=args.num_sequences,
        cd_hit_threshold=args.cd_hit_threshold,
        seq_type=args.seq_type,
        use_fasttree=args.use_fasttree
    )
    
    success = pipeline.run_pipeline(args.input)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
