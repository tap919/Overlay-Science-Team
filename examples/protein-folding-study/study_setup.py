#!/usr/bin/env python3
"""
Protein Folding Study Setup
Initialize study workspace and configuration
"""

from pathlib import Path
import json
from datetime import datetime

# Study configuration
STUDY_ID = 'protein-folding-2026'
STUDY_NAME = 'Protein Folding Analysis 2026'

def setup_study():
    """Initialize study directories and metadata"""
    
    # Create workspace structure
    workspace = Path('workspace')
    input_dir = workspace / 'input' / STUDY_ID
    output_dir = workspace / 'deliverables' / STUDY_ID
    logs_dir = workspace / 'logs'
    
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    
    # Study metadata
    metadata = {
        'study_id': STUDY_ID,
        'study_name': STUDY_NAME,
        'created_at': datetime.now().isoformat(),
        'description': 'Analysis of protein folding pathways using quantum chaos simulations',
        'parameters': {
            'chaos_coefficient': 0.7,
            'iterations': 10000,
            'temperature': 300,  # Kelvin
            'proteins': ['1ubq', '1crn', '1gb1']
        },
        'sources': []
    }
    
    # Save metadata
    with open(input_dir / 'metadata.json', 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"✓ Study '{STUDY_NAME}' initialized")
    print(f"✓ Input directory: {input_dir}")
    print(f"✓ Output directory: {output_dir}")
    print(f"\nNext steps:")
    print(f"1. Add protein structure files to: {input_dir}")
    print(f"2. Run: python data_collection.py")
    print(f"3. Start pipeline in VSCodium: Ctrl+Shift+P → 'Science Lab: Start Science Pipeline'")

if __name__ == '__main__':
    setup_study()
