#!/usr/bin/env python3
"""
Download ebolavirus glycoprotein sequences from UniProt
Downloads separately for each species to ensure complete coverage
Configurable species mapping via JSON config file
"""

import urllib.request
import urllib.parse
import sys
import json
from pathlib import Path

# Default species mapping
DEFAULT_SPECIES = {
    'Zaire': 'organism_name:Zaire ebolavirus',
    'Sudan': 'organism_name:Sudan ebolavirus',
    'Bundibugyo': 'organism_name:Bundibugyo ebolavirus',
    'Taï Forest': 'organism_name:Taï Forest ebolavirus',
    'Reston': 'organism_name:Reston ebolavirus',
    'Bombali': 'organism_name:Bombali virus'
}

def load_species_config(config_file):
    """
    Load species mapping from JSON config file
    
    Expected format:
    {
      "Zaire": "organism_name:Zaire ebolavirus",
      "Sudan": "organism_name:Sudan ebolavirus",
      ...
    }
    """
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
        print(f"Loaded species config from: {config_file}\n")
        return config
    except Exception as e:
        print(f"Error loading config file: {str(e)}")
        print("Using default species mapping\n")
        return DEFAULT_SPECIES


def download_species_gp(species_name, species_query, output_file):
    """
    Download glycoprotein sequences for a single ebolavirus species
    
    Args:
        species_name: Display name of the species
        species_query: UniProt query term for this species
        output_file: Output FASTA file path
    """
    
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    
    # Query for this specific species + glycoprotein
    query = f"{species_query} AND protein_name:glycoprotein"
    
    params = {
        'query': query,
        'format': 'fasta',
        'size': '500'
    }
    
    query_string = urllib.parse.urlencode(params)
    url = f"{base_url}?{query_string}"
    
    try:
        print(f"  Downloading {species_name}...", end=" ", flush=True)
        
        req = urllib.request.Request(url)
        req.add_header('User-Agent', 'Python-Bioinformatics-Script/1.0')
        
        with urllib.request.urlopen(req, timeout=30) as response:
            fasta_content = response.read().decode('utf-8')
        
        num_sequences = fasta_content.count('>')
        
        if num_sequences == 0:
            print(f"✗ NO SEQUENCES FOUND")
            return 0
        
        # Write to file
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, 'w') as f:
            f.write(fasta_content)
        
        print(f"✓ {num_sequences} sequences")
        return num_sequences
        
    except urllib.error.URLError as e:
        print(f"✗ ERROR: {str(e)}")
        return 0
    except Exception as e:
        print(f"✗ ERROR: {str(e)}")
        return 0


def download_all_species(output_dir="sequences", species_map=None):
    """
    Download glycoproteins for all species in the mapping
    Saves each to a separate file and creates a combined file
    """
    
    if species_map is None:
        species_map = DEFAULT_SPECIES
    
    print("="*70)
    print("Downloading Ebolavirus Glycoproteins from UniProt")
    print(f"Species: {', '.join(species_map.keys())}")
    print("="*70 + "\n")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    species_files = {}
    total_sequences = 0
    
    # Download each species separately
    for species_abbrev, species_query in species_map.items():
        safe_name = species_abbrev.lower().replace(' ', '_').replace('ï', 'i')
        output_file = output_dir / f"ebolavirus_gp_{safe_name}.fasta"
        num_seqs = download_species_gp(species_abbrev, species_query, str(output_file))
        
        if num_seqs > 0:
            species_files[species_abbrev] = {
                'file': output_file,
                'count': num_seqs
            }
            total_sequences += num_seqs
    
    # Create combined file
    print(f"\n  Creating combined file...", end=" ", flush=True)
    combined_file = output_dir / "ebolavirus_gp_combined.fasta"
    
    with open(combined_file, 'w') as combined:
        for species_abbrev in species_map.keys():
            if species_abbrev in species_files:
                with open(species_files[species_abbrev]['file'], 'r') as f:
                    combined.write(f.read())
    
    print("✓")
    
    # Print summary
    print("\n" + "="*70)
    print("DOWNLOAD SUMMARY")
    print("="*70)
    print(f"\n{'Species':<20} {'Sequences':<15} {'File'}")
    print("-"*70)
    
    for species_abbrev in species_map.keys():
        if species_abbrev in species_files:
            count = species_files[species_abbrev]['count']
            file_path = species_files[species_abbrev]['file'].name
            print(f"{species_abbrev:<20} {count:<15} {file_path}")
        else:
            print(f"{species_abbrev:<20} {'0':<15} ✗ NOT FOUND")
    
    print("-"*70)
    print(f"{'TOTAL':<20} {total_sequences:<15}")
    print(f"\nCombined file: {combined_file}")
    print("="*70)
    
    return total_sequences > 0


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Download ebolavirus glycoprotein sequences from UniProt",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use default species mapping
  python download_gp_uniprot.py -o sequences
  
  # Use custom species config file
  python download_gp_uniprot.py -o sequences -c species_config.json
  
  # Inside Docker
  docker compose run ebolavirus-analysis python scripts/download_gp_uniprot.py \\
    -o sequences

Species Config File Format (JSON):
  {
    "Zaire": "organism_name:Zaire ebolavirus",
    "Sudan": "organism_name:Sudan ebolavirus",
    "Bundibugyo": "organism_name:Bundibugyo ebolavirus",
    ...
  }
        """
    )
    
    parser.add_argument(
        "-o", "--output",
        default="sequences",
        help="Output directory (default: sequences)"
    )
    parser.add_argument(
        "-c", "--config",
        help="JSON config file with species mapping (optional, uses defaults if not provided)"
    )
    
    args = parser.parse_args()
    
    # Load species mapping
    species_map = None
    if args.config:
        if not Path(args.config).exists():
            print(f"ERROR: Config file not found: {args.config}")
            sys.exit(1)
        species_map = load_species_config(args.config)
    
    success = download_all_species(args.output, species_map)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
