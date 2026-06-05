#!/usr/bin/env python3
"""
Analyze ebolavirus subtype coverage in glycoprotein sequences
Configurable species mapping via JSON config file or command line
"""

import sys
import json
from pathlib import Path
from Bio import SeqIO
from collections import defaultdict

# Default species mapping (can be overridden by config file)
DEFAULT_SPECIES = {
    'Zaire': ['Zaire ebolavirus', 'EBOV'],
    'Sudan': ['Sudan ebolavirus', 'SUDV'],
    'Bundibugyo': ['Bundibugyo ebolavirus', 'BDBV'],
    'Taï Forest': ['Taï Forest ebolavirus', 'TAFV'],
    'Reston': ['Reston ebolavirus', 'RESTV'],
    'Bombali': ['Bombali ebolavirus', 'BOMV']
}

def load_species_config(config_file):
    """
    Load species mapping from JSON config file
    
    Expected format:
    {
      "Zaire": ["Zaire ebolavirus", "EBOV"],
      "Sudan": ["Sudan ebolavirus", "SUDV"],
      ...
    }
    """
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
        print(f"Loaded species config from: {config_file}")
        return config
    except Exception as e:
        print(f"Error loading config file: {str(e)}")
        print("Using default species mapping")
        return DEFAULT_SPECIES


def analyze_subtype_coverage(fasta_file, species_map=None):
    """
    Parse FASTA headers and extract ebolavirus subtypes
    Count sequences per subtype using provided species mapping
    """
    
    if species_map is None:
        species_map = DEFAULT_SPECIES
    
    subtype_counts = {sp: [] for sp in species_map.keys()}
    unclassified = []
    
    print("\n" + "="*70)
    print("EBOLAVIRUS SUBTYPE COVERAGE ANALYSIS")
    print("="*70 + "\n")
    
    # Parse sequences and classify by species
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description
        found = False
        
        # Loop through each species and its keywords
        for species_name, keywords in species_map.items():
            for keyword in keywords:
                if keyword.lower() in header.lower():
                    subtype_counts[species_name].append({
                        'id': record.id,
                        'description': header,
                        'length': len(record.seq)
                    })
                    found = True
                    break
            if found:
                break
        
        if not found:
            unclassified.append({
                'id': record.id,
                'description': header,
                'length': len(record.seq)
            })
    
    # Print summary table
    print("SUBTYPE BREAKDOWN:")
    print("-" * 70)
    
    total_seqs = 0
    for species_name in species_map.keys():
        count = len(subtype_counts[species_name])
        total_seqs += count
        status = "✓" if count > 0 else "✗"
        print(f"{status} {species_name:15} : {count:3} sequences")
    
    print("-" * 70)
    print(f"{'Total classified':15} : {total_seqs:3} sequences")
    print(f"{'Unclassified':15} : {len(unclassified):3} sequences")
    print(f"{'Total':15} : {total_seqs + len(unclassified):3} sequences\n")
    
    # Detailed output by species
    print("\nDETAILED SEQUENCE LIST BY SPECIES:")
    print("=" * 70)
    
    for species_name in species_map.keys():
        seqs = subtype_counts[species_name]
        if seqs:
            print(f"\n{species_name.upper()} EBOLAVIRUS ({len(seqs)} sequences)")
            print("-" * 70)
            for seq_info in seqs:
                desc = seq_info['description']
                print(f"  {seq_info['id']:20} | {seq_info['length']:5} aa | {desc[:60]}")
    
    # Show unclassified sequences if any
    if unclassified:
        print(f"\nUNCLASSIFIED SEQUENCES ({len(unclassified)} sequences)")
        print("-" * 70)
        for seq_info in unclassified:
            desc = seq_info['description']
            print(f"  {seq_info['id']:20} | {seq_info['length']:5} aa | {desc[:60]}")
    
    # Highlight missing species
    missing = [sp for sp in species_map.keys() if len(subtype_counts[sp]) == 0]
    
    if missing:
        print("\n" + "="*70)
        print("WARNING: Missing sequences for the following species:")
        print("="*70)
        for sp in missing:
            print(f"  ✗ {sp}")
        print("\nYou may want to manually add sequences for these species")
        print("from NCBI GenBank or UniProt")
    
    print("\n" + "="*70)
    
    return {
        'subtype_counts': subtype_counts,
        'unclassified_count': len(unclassified),
        'total_sequences': total_seqs + len(unclassified),
        'missing_species': missing
    }


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Analyze ebolavirus subtype coverage in sequence file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use default species mapping
  python analyze_subtypes.py sequences
  
  # Use custom species config file
  python analyze_subtypes.py sequences -c species_config.json
  
  # Inside Docker
  docker compose run ebolavirus-analysis python scripts/analyze_subtypes.py \\
    sequences/ebolavirus_gp_combined.fasta

Species Config File Format (JSON):
  {
    "Zaire": ["Zaire ebolavirus", "EBOV"],
    "Sudan": ["Sudan ebolavirus", "SUDV"],
    "Bundibugyo": ["Bundibugyo ebolavirus", "BDBV"],
    "Taï Forest": ["Taï Forest ebolavirus", "TAFV"],
    "Reston": ["Reston ebolavirus", "RESTV"],
    "Bombali": ["Bombali ebolavirus", "BOMV"]
  }
        """
    )
    
    parser.add_argument("fasta_file", help="Input FASTA file")
    parser.add_argument(
        "-c", "--config",
        help="JSON config file with species mapping (optional, uses defaults if not provided)"
    )
    
    args = parser.parse_args()
    
    if not Path(args.fasta_file).exists():
        print(f"ERROR: File not found: {args.fasta_file}")
        sys.exit(1)
    
    # Load species mapping
    species_map = None
    if args.config:
        if not Path(args.config).exists():
            print(f"ERROR: Config file not found: {args.config}")
            sys.exit(1)
        species_map = load_species_config(args.config)
    
    stats = analyze_subtype_coverage(args.fasta_file, species_map)
    
    # Return non-zero exit if any species is missing
    sys.exit(1 if stats['missing_species'] else 0)


if __name__ == "__main__":
    main()
