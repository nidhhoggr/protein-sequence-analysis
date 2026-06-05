#!/usr/bin/env python3
"""
Post-pipeline analysis utilities for ebolavirus sequences
Visualization, statistics, and metadata extraction
"""

import argparse
import pandas as pd
from pathlib import Path
from Bio import SeqIO, Phylo, AlignIO
from collections import Counter
import json


def analyze_alignment(alignment_file):
    """Analyze multiple sequence alignment statistics"""
    print(f"\n{'='*60}")
    print("Alignment Analysis")
    print(f"{'='*60}")
    
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        
        print(f"Number of sequences: {len(alignment)}")
        print(f"Alignment length: {alignment.get_alignment_length()} bp")
        
        # Check for gaps
        total_gaps = 0
        for seq in alignment:
            total_gaps += str(seq.seq).count("-")
        
        total_length = len(alignment) * alignment.get_alignment_length()
        gap_percentage = (total_gaps / total_length * 100) if total_length > 0 else 0
        
        print(f"Total gaps: {total_gaps} ({gap_percentage:.1f}%)")
        
        # Conservation analysis
        conserved = 0
        for position in range(alignment.get_alignment_length()):
            column = [alignment[i, position] for i in range(len(alignment))]
            column_counts = Counter(column)
            if column_counts.most_common(1)[0][1] >= len(alignment) * 0.9:
                conserved += 1
        
        print(f"Highly conserved positions (≥90%): {conserved}")
        
    except Exception as e:
        print(f"ERROR: {str(e)}")


def analyze_tree(tree_file):
    """Analyze phylogenetic tree structure"""
    print(f"\n{'='*60}")
    print("Tree Analysis")
    print(f"{'='*60}")
    
    try:
        tree = Phylo.read(tree_file, "newick")
        
        # Count terminals
        terminals = [c for c in tree.find_clades(terminal=True)]
        print(f"Number of sequences (terminals): {len(terminals)}")
        
        # Count internal nodes
        internals = [c for c in tree.find_clades(terminal=False)]
        print(f"Number of internal nodes: {len(internals)}")
        
        # Tree depth
        print(f"Tree depth: {tree.distance(terminals[0]):.4f}")
        
        # Clade balance
        def count_leaves(clade):
            if clade.is_terminal():
                return 1
            return sum(count_leaves(c) for c in clade.clades)
        
        clade_sizes = []
        for clade in tree.find_clades(terminal=False):
            clade_sizes.append(count_leaves(clade))
        
        if clade_sizes:
            print(f"Average clade size: {sum(clade_sizes)/len(clade_sizes):.1f} sequences")
        
        print("\nTree structure:")
        Phylo.draw_ascii(tree)
        
    except Exception as e:
        print(f"ERROR: {str(e)}")


def extract_sequences_by_name(fasta_file, names_list):
    """Extract specific sequences by name"""
    print(f"\nExtracting {len(names_list)} sequences...")
    
    output_file = Path(fasta_file).parent / f"extracted_{Path(fasta_file).name}"
    count = 0
    
    with open(output_file, "w") as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in names_list or record.description in names_list:
                SeqIO.write(record, out, "fasta")
                count += 1
    
    print(f"Extracted {count} sequences to {output_file}")
    return output_file


def sequence_statistics(fasta_file):
    """Generate sequence statistics"""
    print(f"\n{'='*60}")
    print("Sequence Statistics")
    print(f"{'='*60}")
    
    lengths = []
    gc_contents = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        lengths.append(len(seq))
        
        # GC content
        gc = (seq.count('G') + seq.count('C')) / len(seq) * 100 if len(seq) > 0 else 0
        gc_contents.append(gc)
    
    if lengths:
        print(f"Sequence count: {len(lengths)}")
        print(f"Length range: {min(lengths)} - {max(lengths)} bp")
        print(f"Mean length: {sum(lengths)/len(lengths):.0f} bp")
        print(f"GC content range: {min(gc_contents):.1f}% - {max(gc_contents):.1f}%")
        print(f"Mean GC content: {sum(gc_contents)/len(gc_contents):.1f}%")


def generate_report(results_dir):
    """Generate comprehensive analysis report"""
    results_dir = Path(results_dir)
    report = []
    
    report.append("=" * 70)
    report.append("EBOLAVIRUS SEQUENCE DIVERSITY ANALYSIS REPORT")
    report.append("=" * 70)
    
    # Check for results files
    files_found = {
        'input': results_dir / "sequences_nonredundant.fasta",
        'aligned': results_dir / "sequences_aligned.fasta",
        'tree': results_dir / "tree.treefile",
        'dataset': results_dir / "ebolavirus_diverse_dataset.fasta",
    }
    
    report.append("\nFiles Generated:")
    for name, path in files_found.items():
        exists = "✓" if path.exists() else "✗"
        report.append(f"  {exists} {name}: {path.name}")
    
    # Sequence statistics
    if files_found['dataset'].exists():
        report.append("\nFinal Dataset Statistics:")
        for record in SeqIO.parse(files_found['dataset'], "fasta"):
            report.append(f"  - {len(list(SeqIO.parse(files_found['dataset'], 'fasta')))} sequences")
            break
    
    report_text = "\n".join(report)
    print(report_text)
    
    # Save report
    report_file = results_dir / "analysis_report.txt"
    with open(report_file, "w") as f:
        f.write(report_text)
    
    print(f"\nReport saved to: {report_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Ebolavirus Sequence Analysis Utilities",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze alignment
  python utils.py --alignment results/sequences_aligned.fasta
  
  # Analyze tree
  python utils.py --tree results/tree.treefile
  
  # Sequence statistics
  python utils.py --stats results/ebolavirus_diverse_dataset.fasta
  
  # Generate full report
  python utils.py --report results/
        """
    )
    
    parser.add_argument("--alignment", help="Analyze multiple sequence alignment")
    parser.add_argument("--tree", help="Analyze phylogenetic tree")
    parser.add_argument("--stats", help="Generate sequence statistics")
    parser.add_argument("--report", help="Generate comprehensive report for results directory")
    
    args = parser.parse_args()
    
    if args.alignment:
        analyze_alignment(args.alignment)
    elif args.tree:
        analyze_tree(args.tree)
    elif args.stats:
        sequence_statistics(args.stats)
    elif args.report:
        generate_report(args.report)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
