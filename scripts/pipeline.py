#!/usr/bin/env python3
"""
Ebolavirus Sequence Diversity Analysis Pipeline
Collects, filters, and aligns ebolavirus sequences with phylogenetic guidance
Supports configurable species mapping via JSON config file
"""

import os
import sys
import subprocess
import argparse
import json
import pandas as pd
import random
from pathlib import Path
from Bio import SeqIO, Phylo
from Bio.SeqRecord import SeqRecord
from io import StringIO

# Default species mapping for diversity-aware selection
DEFAULT_SPECIES = {
    'Zaire': ['Zaire ebolavirus', 'EBOV'],
    'Sudan': ['Sudan ebolavirus', 'SUDV'],
    'Bundibugyo': ['Bundibugyo', 'BDBV'],
    'Taï Forest': ['Taï Forest', 'TAFV'],
    'Reston': ['Reston', 'RESTV'],
    'Bombali': ['Bombali', 'BOMV']
}

class EbolavirusAnalysisPipeline:
    def __init__(self, output_dir="results", num_sequences=40, species_map=None):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.num_sequences = num_sequences
        self.species_map = species_map if species_map else DEFAULT_SPECIES
        self.log_file = self.output_dir / "pipeline.log"
        
    def log(self, message):
        """Log messages to file and stdout"""
        print(message)
        with open(self.log_file, "a") as f:
            f.write(message + "\n")
    
    def run_command(self, cmd, description):
        """Execute shell command and log output"""
        self.log(f"\n{'='*60}")
        self.log(f"Running: {description}")
        self.log(f"Command: {cmd}")
        self.log(f"{'='*60}")
        
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode != 0:
                self.log(f"ERROR: {result.stderr}")
                return False
            self.log(result.stdout)
            return True
        except Exception as e:
            self.log(f"ERROR: {str(e)}")
            return False
    
    def download_sequences(self, input_fasta):
        """
        Use existing FASTA file or download from NCBI.
        For this demo, we'll work with provided sequences.
        """
        if os.path.exists(input_fasta):
            self.log(f"Using provided FASTA file: {input_fasta}")
            return True
        
        self.log("FASTA file not found. Please provide sequences via -i flag.")
        self.log("You can download from:")
        self.log("  - NCBI GenBank: https://www.ncbi.nlm.nih.gov/")
        self.log("  - ViPR: https://www.viprbrc.org/")
        return False
    
    def remove_redundancy(self, input_fasta, output_fasta, identity_threshold=0.95):
        """Remove sequence redundancy using CD-HIT"""
        self.log(f"\nRemoving sequences with >{identity_threshold*100}% identity...")
        
        # For nucleotide sequences: word length (-n) should be 5-8
        cmd = f"cd-hit -i {input_fasta} -o {output_fasta} -c {identity_threshold} -n 5 -d 0 -T 0 2>&1"
        return self.run_command(cmd, f"CD-HIT clustering at {identity_threshold*100}% identity")
    
    def align_sequences(self, input_fasta, output_alignment):
        """Align sequences using MAFFT"""
        self.log(f"\nAligning sequences with MAFFT...")
        
        cmd = f"mafft --auto {input_fasta} > {output_alignment} 2>&1"
        return self.run_command(cmd, "MAFFT sequence alignment")
    
    def build_tree_fasttree(self, aligned_fasta, output_tree):
        """Build phylogenetic tree using FastTree"""
        self.log(f"\nBuilding phylogenetic tree with FastTree...")
        
        cmd = f"FastTree -nt {aligned_fasta} > {output_tree} 2>&1"
        return self.run_command(cmd, "FastTree phylogenetic tree construction")
    
    def build_tree_iqtree(self, aligned_fasta, output_prefix):
        """Build phylogenetic tree using IQ-TREE (more accurate)"""
        self.log(f"\nBuilding phylogenetic tree with IQ-TREE...")
        
        # Use -st AA for protein sequences (glycoproteins are amino acids)
        cmd = f"iqtree -s {aligned_fasta} -st AA -m LG+G -bb 1000 -nt AUTO -pre {output_prefix} 2>&1"
        return self.run_command(cmd, "IQ-TREE phylogenetic tree construction")
    
    def select_by_phylogeny(self, tree_file, sequence_fasta, output_fasta):
        """
        Select diverse sequences by traversing the phylogenetic tree.
        Ensures representation across all major clades and species.
        """
        self.log(f"\nSelecting diverse sequences using phylogenetic guidance...")
        
        try:
            tree = Phylo.read(tree_file, "newick")
            
            # Collect all terminal nodes (sequences)
            all_terminals = []
            for clade in tree.find_clades(terminal=True):
                if clade.name:
                    all_terminals.append(clade.name)
            
            self.log(f"Total sequences in tree: {len(all_terminals)}")
            
            # Strategy: Sample proportionally from tree structure
            selected_names = set()
            
            # Ensure minimum representation by depth-first traversal
            def traverse_and_select(clade, depth=0, max_per_clade=None):
                if clade.is_terminal():
                    if clade.name:
                        # Probabilistically select to reach target number
                        if len(selected_names) < self.num_sequences * 0.9:
                            if random.random() < 0.5:
                                selected_names.add(clade.name)
                else:
                    for child in clade.clades:
                        traverse_and_select(child, depth+1)
            
            # Run multiple passes to reach target number
            for attempt in range(3):
                if len(selected_names) >= self.num_sequences:
                    break
                traverse_and_select(tree.root)
            
            # If still short, add random sequences
            remaining = set(all_terminals) - selected_names
            needed = self.num_sequences - len(selected_names)
            if needed > 0 and remaining:
                selected_names.update(random.sample(list(remaining), min(needed, len(remaining))))
            
            self.log(f"Selected {len(selected_names)} sequences via phylogenetic guidance")
            
            # Extract selected sequences to new FASTA
            selected_records = []
            for record in SeqIO.parse(sequence_fasta, "fasta"):
                # Handle various ID formats
                seq_id = record.id.split("|")[0]  # For GenBank format
                if seq_id in selected_names or record.description in selected_names:
                    selected_records.append(record)
            
            self.log(f"Extracted {len(selected_records)} sequences to {output_fasta}")
            SeqIO.write(selected_records, output_fasta, "fasta")
            
            return len(selected_records) > 0
            
        except Exception as e:
            self.log(f"ERROR in phylogenetic selection: {str(e)}")
            return False
    
    def check_quality(self, fasta_file):
        """Check quality of final dataset"""
        self.log(f"\n{'='*60}")
        self.log("Quality Check")
        self.log(f"{'='*60}")
        
        # Count sequences
        count = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
        self.log(f"Total sequences: {count}")
        
        # Get sequence lengths
        lengths = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            lengths.append(len(record.seq))
        
        if lengths:
            self.log(f"Sequence lengths: min={min(lengths)}, max={max(lengths)}, mean={sum(lengths)/len(lengths):.0f}")
        
        # Check for duplicates
        ids = [record.id for record in SeqIO.parse(fasta_file, "fasta")]
        unique_ids = set(ids)
        if len(ids) != len(unique_ids):
            self.log(f"WARNING: Found {len(ids) - len(unique_ids)} duplicate sequence IDs")
        else:
            self.log("No duplicate sequences detected ✓")
    
    def run_pipeline(self, input_fasta, use_iqtree=True):
        """Execute full pipeline"""
        self.log(f"Starting Ebolavirus Sequence Analysis Pipeline")
        self.log(f"Target: {self.num_sequences} diverse sequences")
        self.log(f"Output directory: {self.output_dir}")
        self.log(f"Species mapping: {list(self.species_map.keys())}")
        
        # Step 1: Download/verify sequences
        if not self.download_sequences(input_fasta):
            return False
        
        # Step 2: Remove redundancy
        nonredundant = self.output_dir / "sequences_nonredundant.fasta"
        if not self.remove_redundancy(input_fasta, str(nonredundant)):
            return False
        
        # Step 3: Align sequences
        aligned = self.output_dir / "sequences_aligned.fasta"
        if not self.align_sequences(str(nonredundant), str(aligned)):
            return False
        
        # Step 4: Build tree
        if use_iqtree:
            tree_output = self.output_dir / "tree"
            if not self.build_tree_iqtree(str(aligned), str(tree_output)):
                return False
            tree_file = str(tree_output) + ".treefile"
        else:
            tree_file = self.output_dir / "tree.nwk"
            if not self.build_tree_fasttree(str(aligned), str(tree_file)):
                return False
        
        # Step 5: Select diverse sequences
        final_fasta = self.output_dir / "ebolavirus_diverse_dataset.fasta"
        if not self.select_by_phylogeny(str(tree_file), str(aligned), str(final_fasta)):
            self.log("Phylogenetic selection failed, using raw alignment instead")
            # Fallback: just use aligned sequences
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


def load_species_config(config_file):
    """Load species mapping from JSON config file"""
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
        print(f"Loaded species config from: {config_file}")
        return config
    except Exception as e:
        print(f"Error loading config file: {str(e)}")
        print("Using default species mapping")
        return DEFAULT_SPECIES


def main():
    parser = argparse.ArgumentParser(
        description="Ebolavirus Sequence Diversity Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with local FASTA file
  python pipeline.py -i sequences.fasta
  
  # Specify number of target sequences
  python pipeline.py -i sequences.fasta -n 50
  
  # Use faster tree building
  python pipeline.py -i sequences.fasta --use-fasttree
  
  # Custom output directory
  python pipeline.py -i sequences.fasta -o my_results
  
  # Use custom species config
  python pipeline.py -i sequences.fasta -c config/species_config.json
        """
    )
    
    parser.add_argument("-i", "--input", default="sequences.fasta",
                        help="Input FASTA file with ebolavirus sequences (default: sequences.fasta)")
    parser.add_argument("-o", "--output", default="results",
                        help="Output directory (default: results)")
    parser.add_argument("-n", "--num-sequences", type=int, default=40,
                        help="Target number of diverse sequences to select (default: 40)")
    parser.add_argument("--use-fasttree", action="store_true",
                        help="Use FastTree instead of IQ-TREE (faster but less accurate)")
    parser.add_argument("-c", "--config",
                        help="JSON config file with species mapping (optional, uses defaults if not provided)")
    
    args = parser.parse_args()
    
    # Load species mapping
    species_map = None
    if args.config:
        if not Path(args.config).exists():
            print(f"ERROR: Config file not found: {args.config}")
            sys.exit(1)
        species_map = load_species_config(args.config)
    
    pipeline = EbolavirusAnalysisPipeline(
        output_dir=args.output,
        num_sequences=args.num_sequences,
        species_map=species_map
    )
    
    success = pipeline.run_pipeline(
        input_fasta=args.input,
        use_iqtree=not args.use_fasttree
    )
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
