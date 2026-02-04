# Overlay Science Lab - VSCodium Extension

A powerful VSCodium extension for the Overlay Science Team that brings enhanced scientific research capabilities directly into your IDE.

## Features

### 🔬 Integrated Science Pipeline
- **One-Click Pipeline Execution**: Start multi-agent scientific workflows directly from VSCodium
- **Real-Time Progress Tracking**: Monitor pipeline stages with live updates in the status bar
- **Study Management**: Create, organize, and track multiple research studies

### 🧬 Scientific Tools Integration
- **GLCCE Quantum Chaos**: Run quantum simulations with configurable parameters
- **Digital Lab**: Execute biotech simulations for protein studies
- **Book Lab**: Generate research papers and book chapters
- **Comic Engine**: Create visual metaphors for complex concepts
- **DataLite**: Perform analytics and visualization

### 📊 Enhanced Research Capabilities
- **Scientific API Snippets**: Quick access to 7+ popular scientific APIs
  - PubChemPy for compound data
  - Biopython for PubMed queries
  - BioThings for gene information
  - PyBioPortal for cancer genomics
  - ChEMBL for drug discovery
  - Scikit-Bio for sequence analysis
  - BioMart for gene attributes

### 📝 Science Notebooks
- **Interactive Notebooks**: Create `.science-nb` notebooks for exploratory research
- **Python Integration**: Full Python support with scientific libraries
- **Live Code Execution**: Run code cells and see results instantly

### 👁️ Visual Dashboard
- **Pipeline Visualization**: See all 6 pipeline stages at a glance
- **Deliverables View**: Quick access to generated papers, chapters, and analytics
- **Agent Status**: Monitor the health and performance of all pipeline agents

## Installation

### From Source
1. Clone this repository
2. Open in VSCodium
3. Press F5 to run in development mode

### From VSIX (when available)
1. Download the `.vsix` file
2. Open VSCodium
3. Go to Extensions → Install from VSIX
4. Select the downloaded file

## Quick Start

### 1. Configure Your Environment

Open VSCodium Settings (Ctrl+,) and search for "Overlay Science Lab":

```json
{
  "overlay-science-lab.workspacePath": "./science_workspace",
  "overlay-science-lab.apiKey": "your-anthropic-api-key",
  "overlay-science-lab.glccePath": "/path/to/glcce",
  "overlay-science-lab.digitalLabPath": "/path/to/digital-lab",
  "overlay-science-lab.autoStartPipeline": false,
  "overlay-science-lab.enableRealTimeUpdates": true
}
```

### 2. Create Your First Study

**Option A: Command Palette**
1. Press `Ctrl+Shift+P`
2. Type "Science Lab: Create New Study"
3. Enter study name
4. Add source files to the created directory

**Option B: Sidebar**
1. Click the Science Lab icon in the Activity Bar
2. Click the "+" button in Active Studies view
3. Enter study name
4. Upload your research sources

### 3. Run the Pipeline

**Option A: Command Palette**
- Press `Ctrl+Shift+P`
- Type "Science Lab: Start Science Pipeline"
- Enter study ID

**Option B: Status Bar**
- Click the "Science Lab: Ready" item in the status bar
- Select your study from the dashboard

### 4. View Results

Generated deliverables appear in:
- `science_workspace/deliverables/<study_id>/`

Access them via:
- Deliverables tree view in the sidebar
- "Science Lab: View Deliverables" command

## Commands

All commands are prefixed with "Science Lab:" in the Command Palette:

| Command | Description | Shortcut |
|---------|-------------|----------|
| Start Science Pipeline | Execute full pipeline for a study | - |
| Open Science Dashboard | View pipeline status and metrics | - |
| Create New Study | Initialize a new research study | - |
| Upload Research Sources | Add files to a study | - |
| View Deliverables | Open deliverables directory | - |
| Run Quantum Chaos Simulation | Execute GLCCE simulation | - |
| Run Biotech Simulation | Execute Digital Lab analysis | - |
| Open Science Notebook | Create interactive notebook | - |
| Insert Scientific API Snippet | Add pre-built API code | - |
| Stop Active Pipeline | Cancel running pipeline | - |

## Code Snippets

Type these prefixes in Python files to insert scientific code:

### API Snippets
- `pubchem-search` - Search compounds in PubChem
- `biopython-pubmed` - Query PubMed articles
- `biothings-gene` - Fetch gene information
- `science-pipeline-setup` - Initialize study structure

### Configuration Snippets
- `quantum-sim-params` - Quantum simulation parameters
- `data-analysis-pipeline` - Complete data analysis workflow

## Extension Views

### Active Studies
- Lists all running and completed studies
- Shows progress for each pipeline stage
- Quick actions: Open, Stop, View Results

### Deliverables
- Research papers (PDF/LaTeX)
- Book chapters (PDF/Markdown)
- Analytics dashboards (HTML)
- Metaphor images (PNG)

### Scientific Tools
- Quick access to all integrated tools
- Tool status indicators
- Direct launch capabilities

### Pipeline Agents
- 6 specialized agents:
  - Source Ingester
  - Quantum Runner
  - Biotech Simulator
  - Analytics Compiler
  - Writer Agent
  - Metaphor Agent
- Real-time status updates
- Performance metrics

## Pipeline Stages

### Stage 1: Source Ingestion (📥)
- Parses PDFs, CSVs, and research materials
- Extracts key information and hypotheses
- Prepares data for analysis

### Stage 2: Quantum Chaos Analysis (⚛️)
- Runs GLCCE simulations
- Analyzes chaos coefficients
- Generates quantum insights

### Stage 3: Biotech Simulation (🧬)
- Digital Lab protein studies
- Pathway analysis
- Molecular dynamics

### Stage 4: Insight Synthesis (📊)
- DataLite analytics
- Statistical analysis
- Visualization generation

### Stage 5: Paper & Book Writing (📝)
- Generates research papers
- Creates book chapters
- Formats citations

### Stage 6: Metaphor Enhancement (🎨)
- Creates visual metaphors
- Enhances understanding
- Generates illustrations

## Configuration Options

### Workspace Settings

| Setting | Type | Default | Description |
|---------|------|---------|-------------|
| `workspacePath` | string | `./science_workspace` | Path to science workspace |
| `apiKey` | string | `""` | Anthropic API key |
| `glccePath` | string | `""` | GLCCE executable path |
| `digitalLabPath` | string | `""` | Digital Lab path |
| `bookLabPath` | string | `""` | Book Lab path |
| `comicEnginePath` | string | `""` | Comic Engine path |
| `autoStartPipeline` | boolean | `false` | Auto-start on file upload |
| `enableRealTimeUpdates` | boolean | `true` | Show live progress |
| `notebookIntegration` | boolean | `true` | Enable notebooks |
| `showAgentMetrics` | boolean | `true` | Display agent metrics |

## Troubleshooting

### Extension Not Activating
- Check that you have a workspace folder open
- Verify VSCodium version >= 1.80.0
- Check the Output panel (View → Output → Overlay Science Lab)

### Pipeline Fails to Start
- Verify API key is configured
- Check tool paths in settings
- Ensure `science_workspace` directory exists
- Review logs in `science_workspace/logs/`

### API Snippets Not Working
- Verify you're in a Python file
- Type the full prefix and wait for autocomplete
- Check that Python extension is installed

### Status Bar Not Updating
- Enable "Real-Time Updates" in settings
- Restart VSCodium
- Check for conflicting extensions

## Development

### Building from Source

```bash
cd vscodium-extension
npm install
npm run lint
```

### Running Tests

```bash
npm test
```

### Debugging

1. Open the extension folder in VSCodium
2. Press F5 to launch Extension Development Host
3. Set breakpoints in `src/extension.js`

## Requirements

- VSCodium / VSCode >= 1.80.0
- Python 3.8+ (for scientific APIs)
- Node.js 18+ (for extension development)

### Optional Dependencies
- Anthropic API key (for AI agents)
- GLCCE (for quantum simulations)
- Digital Lab (for biotech studies)
- Book Lab (for paper generation)
- Comic Engine (for metaphors)

## Roadmap

- [ ] Real-time collaboration features
- [ ] Voice control integration
- [ ] Mobile companion app
- [ ] Auto-submission to preprint servers
- [ ] Peer review simulation
- [ ] Citation network analysis
- [ ] GPU acceleration support
- [ ] Cloud pipeline execution

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## License

MIT License - See LICENSE file for details

## Support

- **Issues**: [GitHub Issues](https://github.com/tap919/Overlay-Science-Team/issues)
- **Discussions**: [GitHub Discussions](https://github.com/tap919/Overlay-Science-Team/discussions)
- **Documentation**: [Wiki](https://github.com/tap919/Overlay-Science-Team/wiki)

## Acknowledgments

Built with:
- VSCodium Extension API
- Node.js
- Anthropic Claude API
- Scientific Python ecosystem

---

**Happy Researching! 🔬**
