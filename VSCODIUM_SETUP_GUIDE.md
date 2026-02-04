# VSCodium Science Lab IDE - Installation & Setup Guide

## Overview

This guide will help you set up the complete VSCodium-based Science Lab IDE with all enhanced scientific capabilities.

## Prerequisites

### Required Software
- **VSCodium** 1.80.0 or later (or VSCode)
- **Python** 3.8 or later
- **Node.js** 18.x or later
- **Git** 2.x or later

### Optional Tools
- GLCCE (Quantum Chaos simulation)
- Digital Lab (Biotech simulation)
- Book Lab (Paper generation)
- Comic Engine (Visual metaphors)

## Installation Steps

### 1. Install VSCodium

#### Linux
```bash
# Ubuntu/Debian
wget -qO - https://gitlab.com/paulcarroty/vscodium-deb-rpm-repo/raw/master/pub.gpg \
    | gpg --dearmor \
    | sudo dd of=/usr/share/keyrings/vscodium-archive-keyring.gpg

echo 'deb [ signed-by=/usr/share/keyrings/vscodium-archive-keyring.gpg ] https://download.vscodium.com/debs vscodium main' \
    | sudo tee /etc/apt/sources.list.d/vscodium.list

sudo apt update && sudo apt install codium
```

#### macOS
```bash
brew install --cask vscodium
```

#### Windows
Download from: https://vscodium.com/

### 2. Install Python Scientific Stack

```bash
# Create virtual environment
python3 -m venv ~/science-env
source ~/science-env/bin/activate  # On Windows: science-env\Scripts\activate

# Install core scientific packages
pip install --upgrade pip
pip install numpy pandas matplotlib scipy seaborn jupyter

# Install enhanced scientific APIs
pip install pubchempy biopython biothings-client pybioportal
pip install chembl-downloader rdkit scikit-bio mygene metapub
pip install biopandas prody pysam biomart pyensembl bioservices

# Install data analysis tools
pip install statsmodels scikit-learn plotly dash

# Install molecular dynamics (optional)
pip install MDAnalysis nglview
```

### 3. Install Overlay Science Lab Extension

#### Option A: From Source (Development)
```bash
# Clone the repository
git clone https://github.com/tap919/Overlay-Science-Team.git
cd Overlay-Science-Team/vscodium-extension

# Install dependencies
npm install

# Open in VSCodium and press F5 to launch
codium .
```

#### Option B: From VSIX (Future Release)
1. Download `overlay-science-lab-1.0.0.vsix`
2. Open VSCodium
3. Press Ctrl+Shift+P (Cmd+Shift+P on Mac)
4. Type "Extensions: Install from VSIX"
5. Select the downloaded file

### 4. Configure Extension Settings

Open VSCodium Settings (Ctrl+,) and configure:

```json
{
  // Workspace configuration
  "overlay-science-lab.workspacePath": "./science_workspace",
  
  // API keys
  "overlay-science-lab.apiKey": "sk-ant-your-anthropic-key-here",
  
  // Tool paths (configure if installed)
  "overlay-science-lab.glccePath": "/usr/local/bin/glcce",
  "overlay-science-lab.digitalLabPath": "/opt/digital-lab/bin/digital-lab",
  "overlay-science-lab.bookLabPath": "/opt/book-lab/book-lab",
  "overlay-science-lab.comicEnginePath": "/opt/comic-engine/comic-engine",
  
  // Pipeline settings
  "overlay-science-lab.autoStartPipeline": false,
  "overlay-science-lab.enableRealTimeUpdates": true,
  "overlay-science-lab.notebookIntegration": true,
  "overlay-science-lab.showAgentMetrics": true
}
```

### 5. Initialize Science Workspace

```bash
# Create workspace structure
mkdir -p ~/science_workspace/{input,deliverables,logs,temp}

# Open workspace in VSCodium
codium ~/science_workspace
```

## Enhanced Scientific Features

### 1. Scientific API Integration

The extension provides 15+ scientific API integrations:

#### Compound & Drug APIs
- **PubChemPy**: Compound information from PubChem
- **ChEMBL**: Drug discovery database
- **RDKit**: Cheminformatics toolkit

#### Literature APIs
- **Biopython**: PubMed, Gene, Protein databases
- **MetaPub**: Enhanced PubMed access

#### Gene & Protein APIs
- **BioThings**: MyGene.info, MyVariant.info
- **MyGene**: Gene annotation service
- **PyBioPortal**: Cancer genomics data

#### Sequence & Structure APIs
- **BioPandas**: Molecular structures in DataFrames
- **ProDy**: Protein dynamics analysis
- **scikit-bio**: Sequence analysis
- **pysam**: Genomic file handling

#### Genomics APIs
- **BioMart**: Ensembl BioMart wrapper
- **PyEnsembl**: Ensembl reference genome
- **BioServices**: 25+ biological web services

### 2. Code Snippets

Type these prefixes in Python files:

| Prefix | Description |
|--------|-------------|
| `pubchem-search` | Search PubChem for compounds |
| `biopython-pubmed` | Query PubMed articles |
| `biothings-gene` | Fetch gene information |
| `science-pipeline-setup` | Initialize study structure |
| `quantum-sim-params` | Quantum simulation parameters |
| `data-analysis-pipeline` | Complete data analysis |

### 3. Integrated Tools

#### GLCCE (Quantum Chaos)
```bash
# Install GLCCE (if available)
# Contact your system administrator or visit project site
```

#### Digital Lab (Biotech Simulation)
```bash
# Install GROMACS (open-source molecular dynamics)
sudo apt-get install gromacs  # Ubuntu/Debian
brew install gromacs           # macOS
```

### 4. Science Notebooks

Create interactive notebooks:
1. Press Ctrl+Shift+P
2. Type "Science Lab: Open Science Notebook"
3. Write Python code with scientific libraries
4. Execute cells interactively

## Usage Examples

### Example 1: Quick PubMed Search

```python
# Type: biopython-pubmed
from Bio import Entrez

Entrez.email = 'your@email.com'
handle = Entrez.esearch(db='pubmed', term='CRISPR', retmax=10)
record = Entrez.read(handle)
print(record['IdList'])
```

### Example 2: Gene Information Query

```python
# Type: biothings-gene
from biothings_client import get_client

mg = get_client('gene')
results = mg.querymany(['TP53', 'BRCA1'], scopes='symbol')
for r in results:
    print(f"{r['symbol']}: {r['name']}")
```

### Example 3: Compound Structure Search

```python
# Type: pubchem-search
import pubchempy as pcp

compounds = pcp.get_compounds('caffeine', 'name')
for c in compounds:
    print(f"Formula: {c.molecular_formula}")
    print(f"MW: {c.molecular_weight}")
    print(f"SMILES: {c.isomeric_smiles}")
```

### Example 4: Full Pipeline Study

```bash
# 1. Create new study
Ctrl+Shift+P → "Science Lab: Create New Study"
# Enter: "protein-folding-2026"

# 2. Add source files
# Copy PDFs, CSVs to: science_workspace/input/protein-folding-2026/

# 3. Start pipeline
Ctrl+Shift+P → "Science Lab: Start Science Pipeline"
# Enter: "protein-folding-2026"

# 4. Monitor progress
# Watch status bar and pipeline dashboard

# 5. View results
Ctrl+Shift+P → "Science Lab: View Deliverables"
```

## Advanced Configuration

### Custom Task Definitions

Add to `.vscode/tasks.json`:

```json
{
  "version": "2.0.0",
  "tasks": [
    {
      "label": "Run Science Pipeline",
      "type": "science-pipeline",
      "command": "start",
      "studyId": "${input:studyId}",
      "problemMatcher": ["$science-pipeline"]
    },
    {
      "label": "Run Quantum Simulation",
      "type": "shell",
      "command": "python",
      "args": [
        "-c",
        "from enhanced_scientific_apis import *; run_quantum_sim()"
      ]
    }
  ],
  "inputs": [
    {
      "id": "studyId",
      "type": "promptString",
      "description": "Enter study ID"
    }
  ]
}
```

### Custom Keybindings

Add to `keybindings.json`:

```json
[
  {
    "key": "ctrl+alt+s",
    "command": "overlay-science-lab.startPipeline",
    "when": "editorTextFocus"
  },
  {
    "key": "ctrl+alt+d",
    "command": "overlay-science-lab.openDashboard"
  },
  {
    "key": "ctrl+alt+n",
    "command": "overlay-science-lab.openNotebook"
  }
]
```

## Troubleshooting

### Extension Not Loading
```bash
# Check VSCodium version
codium --version  # Should be >= 1.80.0

# Check extension host logs
# View → Output → Log (Extension Host)
```

### Python APIs Not Working
```bash
# Verify Python environment
python --version
pip list | grep -E "pubchempy|biopython|biothings"

# Reinstall if needed
pip install --force-reinstall pubchempy biopython biothings-client
```

### Pipeline Fails to Start
```bash
# Check workspace structure
ls -la science_workspace/
# Should have: input/, deliverables/, logs/

# Check logs
cat science_workspace/logs/pipeline.log

# Verify API key
# Settings → Overlay Science Lab → API Key
```

### Status Bar Not Updating
1. Restart VSCodium
2. Check "Enable Real-Time Updates" setting
3. Open Output panel: View → Output → Overlay Science Lab

## Performance Optimization

### For Low-Spec Machines

```json
{
  "overlay-science-lab.autoStartPipeline": false,
  "overlay-science-lab.enableRealTimeUpdates": false,
  "overlay-science-lab.showAgentMetrics": false
}
```

### For High-Performance Workstations

```json
{
  "overlay-science-lab.autoStartPipeline": true,
  "overlay-science-lab.enableRealTimeUpdates": true,
  "overlay-science-lab.showAgentMetrics": true
}
```

## Integration with External Tools

### Jupyter Notebooks
```bash
# Install Jupyter
pip install jupyter notebook

# Launch from VSCodium
# Terminal → New Terminal
jupyter notebook
```

### JupyterLab
```bash
# Install JupyterLab
pip install jupyterlab

# Launch
jupyter lab
```

### RStudio (for R users)
```bash
# Install R kernel
pip install rpy2

# Use R in notebooks
```

## Security Best Practices

1. **API Keys**: Store in environment variables, not in code
2. **Credentials**: Use VSCodium Secret Storage
3. **Data**: Keep sensitive data in encrypted folders
4. **Git**: Add `.env` and API keys to `.gitignore`

```bash
# .gitignore
.env
*.env
science_workspace/logs/
science_workspace/temp/
**/*.key
**/*.pem
```

## Next Steps

1. **Tutorial**: Complete the Quick Start tutorial
2. **Examples**: Try the example studies in `examples/`
3. **Documentation**: Read the full docs at `/docs`
4. **Community**: Join discussions on GitHub

## Support

- **GitHub Issues**: https://github.com/tap919/Overlay-Science-Team/issues
- **Discussions**: https://github.com/tap919/Overlay-Science-Team/discussions
- **Email**: support@overlay-science.org

## License

MIT License - See LICENSE file

---

**Ready to revolutionize your scientific workflow!** 🔬🚀
