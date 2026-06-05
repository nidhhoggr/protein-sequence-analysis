#!/usr/bin/env python3
"""
Download ebolavirus sequences from UniProt
Downloads separately for each species to ensure complete coverage
Configurable species mapping and protein filters via JSON config file
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

# Default protein filter
DEFAULT_PROTEIN_FILTER = 'protein_name:glycoprotein'

def load_config(config_file):
    """
    Load configuration from JSON config file
    
    Expected format:
    {
      "species": {
        "Zaire": "organism_name:Zaire ebolavirus",
        ...
      },
      "protein_filter": "protein_name:glycoprotein"
    }
    """
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        # Check if it's the new format with nested species and protein_filter
        if 'species' in config:
            species_map = config['species']
            protein_filter = config.get('protein_filter', DEFAULT_PROTEIN_FILTER)
        else:
            # Old format: just species mapping
            species_map = config
            protein_filter = DEFAULT_PROTEIN_FILTER
        
        print(f"Loaded config from: {config_file}")
        print(f"  Species: {len(species_map)}")
        print(f"  Protein filter: {protein_filter}\n")
        return species_map, protein_filter
    except Exception as e:
        print(f"Error loading config file: {str(e)}")
        print("Using default species mapping and protein filter\n")
        return DEFAULT_SPECIES, DEFAULT_PROTEIN_FILTER


def download_species_gp(species_name, species_query, protein_filter, output_file):
    """
    Download sequences for a single ebolavirus species
    
    Args:
        species_name: Display name of the species
        species_query: UniProt query term for this species
        protein_filter: UniProt query term for protein type
        output_file: Output FASTA file path
    """
    
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    
    # Query for this specific species + protein filter
    query = f"{species_query} AND {protein_filter}"
    
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


def download_all_species(output_dir="sequences", species_map=None, protein_filter=None):
    """
    Download sequences for all species in the mapping
    Saves each to a separate file and creates a combined file
    """
    
    if species_map is None:
        species_map = DEFAULT_SPECIES
    if protein_filter is None:
        protein_filter = DEFAULT_PROTEIN_FILTER
    
    print("="*70)
    print("Downloading Ebolavirus Sequences from UniProt")
    print(f"Species: {', '.join(species_map.keys())}")
    print(f"Protein filter: {protein_filter}")
    print("="*70 + "\n")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    species_files = {}
    total_sequences = 0
    
    # Download each species separately
    for species_abbrev, species_query in species_map.items():
        safe_name = species_abbrev.lower().replace(' ', '_').replace('ï', 'i')
        output_file = output_dir / f"ebolavirus_gp_{safe_name}.fasta"
        num_seqs = download_species_gp(species_abbrev, species_query, protein_filter, str(output_file))
        
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
        description="Download ebolavirus sequences from UniProt",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use default settings (glycoprotein)
  python download_gp_uniprot.py -o sequences
  
  # Use custom config from config/ directory
  python download_gp_uniprot.py -o sequences -c config/species_config.json
  
  # Inside Docker
  docker compose run ebolavirus-analysis python scripts/download_gp_uniprot.py \\
    -o sequences -c config/species_config.json

Config File Location:
  Place config files in: ./config/
  Access in container as: /data/config/

Config File Format (JSON):
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

Available protein filters:
  - protein_name:glycoprotein
  - protein_name:nucleoprotein
  - protein_name:polymerase
  - (any valid UniProt protein_name filter)
        """
    )
    
    parser.add_argument(
        "-o", "--output",
        default="sequences",
        help="Output directory (default: sequences)"
    )
    parser.add_argument(
        "-c", "--config",
        help="JSON config file path (e.g., config/species_config.json)"
    )
    
    args = parser.parse_args()
    
    # Load configuration
    species_map = None
    protein_filter = None
    
    if args.config:
        if not Path(args.config).exists():
            print(f"ERROR: Config file not found: {args.config}")
            print(f"Make sure the config file is in the ./config/ directory")
            sys.exit(1)
        species_map, protein_filter = load_config(args.config)
    
    success = download_all_species(args.output, species_map, protein_filter)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
