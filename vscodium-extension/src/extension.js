const vscode = require('vscode');
const path = require('path');
const fs = require('fs').promises;
const ScientificAPICompletionProvider = require('./intellisense');
const { runAgenticWorkflow } = require('./agenticWorkflow');

let statusBarItem;
let activePipelines = new Map();
let studyTreeProvider;
let deliverableTreeProvider;
let toolsTreeProvider;
let agentTreeProvider;
let intellisenseProvider;

/**
 * Extension activation entry point
 */
function activate(context) {
    console.log('Overlay Science Lab extension is now active!');

    // Initialize status bar
    statusBarItem = vscode.window.createStatusBarItem(vscode.StatusBarAlignment.Left, 100);
    statusBarItem.text = '$(beaker) Science Lab: Ready';
    statusBarItem.command = 'overlay-science-lab.openDashboard';
    statusBarItem.show();
    context.subscriptions.push(statusBarItem);

    // Initialize tree view providers
    studyTreeProvider = new StudyTreeProvider();
    deliverableTreeProvider = new DeliverableTreeProvider();
    toolsTreeProvider = new ToolsTreeProvider();
    agentTreeProvider = new AgentTreeProvider();

    vscode.window.registerTreeDataProvider('activeStudies', studyTreeProvider);
    vscode.window.registerTreeDataProvider('deliverables', deliverableTreeProvider);
    vscode.window.registerTreeDataProvider('tools', toolsTreeProvider);
    vscode.window.registerTreeDataProvider('agents', agentTreeProvider);

    // Register IntelliSense providers
    intellisenseProvider = new ScientificAPICompletionProvider();
    
    context.subscriptions.push(
        vscode.languages.registerCompletionItemProvider(
            { scheme: 'file', language: 'python' },
            intellisenseProvider,
            '.'
        )
    );
    
    context.subscriptions.push(
        vscode.languages.registerHoverProvider(
            { scheme: 'file', language: 'python' },
            intellisenseProvider
        )
    );
    
    context.subscriptions.push(
        vscode.languages.registerSignatureHelpProvider(
            { scheme: 'file', language: 'python' },
            intellisenseProvider,
            '(', ','
        )
    );

    // Register commands
    registerCommands(context);

    // Check for science workspace and initialize
    initializeWorkspace();

    // Setup file watchers for auto-processing
    setupFileWatchers(context);

    vscode.window.showInformationMessage('Overlay Science Lab initialized successfully!');
}

/**
 * Register all extension commands
 */
function registerCommands(context) {
    // Start Pipeline Command
    context.subscriptions.push(
        vscode.commands.registerCommand('overlay-science-lab.startPipeline', async () => {
            const studyId = await vscode.window.showInputBox({
                prompt: 'Enter study ID or name',
                placeHolder: 'e.g., protein-folding-2026',
                value: `study_${Date.now()}`
            });

            if (studyId) {
                await startPipeline(studyId);
            }
        })
    );

    // Open Dashboard Command
    context.subscriptions.push(
        vscode.commands.registerCommand('overlay-science-lab.openDashboard', () => {
            openDashboardPanel(context);
        })
    );

    // New Study Command
    context.subscriptions.push(
        vscode.commands.registerCommand('overlay-science-lab.newStudy', async () => {
            await createNewStudy();
        })
    );

    // Upload Sources Command
    context.subscriptions.push(
        vscode.commands.registerCommand('overlay-science-lab.uploadSources', async () => {
            await uploadSources();
        })
    );

    // View Deliverables Command
    context.subscriptions.push(
        vscode.commands.registerCommand('overlay-science-lab.viewDeliverables', async () => {
            await viewDeliverables();
        })
    );

    // Run Quantum Simulation Command
    context.subscriptions.push(
        vscode.commands.registerCommand('overlay-science-lab.runQuantumSim', async () => {
            await runQuantumSimulation();
        })
    );

    // Run Biotech Simulation Command
    context.subscriptions.push(
        vscode.commands.registerCommand('overlay-science-lab.runBiotechSim', async () => {
            await runBiotechSimulation();
        })
    );

    // Open Notebook Command
    context.subscriptions.push(
        vscode.commands.registerCommand('overlay-science-lab.openNotebook', async () => {
            await openScienceNotebook();
        })
    );

    // Insert API Snippet Command
    context.subscriptions.push(
        vscode.commands.registerCommand('overlay-science-lab.insertApiSnippet', async () => {
            await insertApiSnippet();
        })
    );

    // Stop Pipeline Command
    context.subscriptions.push(
        vscode.commands.registerCommand('overlay-science-lab.stopPipeline', async () => {
            await stopActivePipeline();
        })
    );

    // Agentic automation workflow (end-to-end orchestration)
    context.subscriptions.push(
        vscode.commands.registerCommand('overlay-science-lab.runAgenticWorkflow', async () => {
            await runAgenticWorkflow({
                statusBarItem,
                studyTreeProvider,
                deliverableTreeProvider
            });
        })
    );
}

/**
 * Initialize science workspace
 */
async function initializeWorkspace() {
    const config = vscode.workspace.getConfiguration('overlay-science-lab');
    const workspacePath = config.get('workspacePath', './science_workspace');
    
    const workspaceFolders = vscode.workspace.workspaceFolders;
    if (!workspaceFolders || workspaceFolders.length === 0) {
        return;
    }

    const rootPath = workspaceFolders[0].uri.fsPath;
    const scienceWorkspace = path.join(rootPath, workspacePath);

    try {
        await fs.mkdir(scienceWorkspace, { recursive: true });
        await fs.mkdir(path.join(scienceWorkspace, 'input'), { recursive: true });
        await fs.mkdir(path.join(scienceWorkspace, 'deliverables'), { recursive: true });
        await fs.mkdir(path.join(scienceWorkspace, 'logs'), { recursive: true });
        
        // Refresh tree views
        studyTreeProvider.refresh();
        deliverableTreeProvider.refresh();
    } catch (error) {
        console.error('Error initializing workspace:', error);
    }
}

/**
 * Setup file watchers for auto-processing
 */
function setupFileWatchers(context) {
    const config = vscode.workspace.getConfiguration('overlay-science-lab');
    const autoStart = config.get('autoStartPipeline', false);

    if (!autoStart) {
        return;
    }

    const workspaceFolders = vscode.workspace.workspaceFolders;
    if (!workspaceFolders) {
        return;
    }

    const workspacePath = config.get('workspacePath', 'science_workspace');
    const pattern = new vscode.RelativePattern(
        workspaceFolders[0],
        `**/${workspacePath}/input/**/*`
    );

    const watcher = vscode.workspace.createFileSystemWatcher(pattern);
    
    watcher.onDidCreate(async (uri) => {
        const studyId = `study_${Date.now()}`;
        vscode.window.showInformationMessage(
            `New source detected. Start pipeline for study ${studyId}?`,
            'Yes', 'No'
        ).then(answer => {
            if (answer === 'Yes') {
                startPipeline(studyId);
            }
        });
    });

    context.subscriptions.push(watcher);
}

/**
 * Start the science pipeline
 */
async function startPipeline(studyId) {
    const workspaceFolders = vscode.workspace.workspaceFolders;
    
    if (!workspaceFolders) {
        vscode.window.showErrorMessage('No workspace folder open');
        return;
    }
    
    // Register pipeline as active
    activePipelines.set(studyId, { startTime: Date.now() });
    
    statusBarItem.text = `$(sync~spin) Science Lab: Running ${studyId}`;
    
    vscode.window.withProgress({
        location: vscode.ProgressLocation.Notification,
        title: `Running Science Pipeline: ${studyId}`,
        cancellable: true
    }, async (progress, token) => {
        progress.report({ increment: 0, message: 'Initializing...' });
        
        // Simulate pipeline stages
        const stages = [
            { name: 'Source Ingestion', duration: 2000 },
            { name: 'Quantum Chaos Analysis', duration: 3000 },
            { name: 'Biotech Simulation', duration: 3000 },
            { name: 'Insight Synthesis', duration: 2000 },
            { name: 'Paper Writing', duration: 2000 },
            { name: 'Metaphor Enhancement', duration: 2000 }
        ];

        let completed = 0;
        for (const stage of stages) {
            if (token.isCancellationRequested) {
                break;
            }
            
            progress.report({ 
                increment: 100 / stages.length, 
                message: stage.name 
            });
            
            await new Promise(resolve => setTimeout(resolve, stage.duration));
            completed++;
        }

        if (completed === stages.length) {
            vscode.window.showInformationMessage(
                `Pipeline completed for ${studyId}! Deliverables ready.`
            );
            statusBarItem.text = '$(check) Science Lab: Complete';
            
            // Refresh tree views
            studyTreeProvider.refresh();
            deliverableTreeProvider.refresh();
        } else {
            statusBarItem.text = '$(warning) Science Lab: Cancelled';
        }

        // Remove from active pipelines
        activePipelines.delete(studyId);

        setTimeout(() => {
            statusBarItem.text = '$(beaker) Science Lab: Ready';
        }, 3000);
    });
}

/**
 * Open dashboard webview panel
 */
function openDashboardPanel(context) {
    const panel = vscode.window.createWebviewPanel(
        'scienceDashboard',
        'Science Lab Dashboard',
        vscode.ViewColumn.One,
        {
            enableScripts: true,
            retainContextWhenHidden: true
        }
    );

    panel.webview.html = getDashboardHtml();
}

/**
 * Get HTML content for dashboard
 */
function getDashboardHtml() {
    return `<!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Science Lab Dashboard</title>
        <style>
            body {
                font-family: var(--vscode-font-family);
                background-color: var(--vscode-editor-background);
                color: var(--vscode-editor-foreground);
                padding: 20px;
            }
            .dashboard-header {
                margin-bottom: 30px;
            }
            .dashboard-header h1 {
                margin: 0 0 10px 0;
                color: var(--vscode-textLink-foreground);
            }
            .stats-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 20px;
                margin-bottom: 30px;
            }
            .stat-card {
                background: var(--vscode-editor-inactiveSelectionBackground);
                padding: 20px;
                border-radius: 8px;
                border: 1px solid var(--vscode-panel-border);
            }
            .stat-value {
                font-size: 32px;
                font-weight: bold;
                margin-bottom: 5px;
            }
            .stat-label {
                font-size: 14px;
                opacity: 0.8;
            }
            .pipeline-stages {
                margin-top: 30px;
            }
            .stage {
                padding: 15px;
                margin-bottom: 10px;
                background: var(--vscode-editor-inactiveSelectionBackground);
                border-radius: 8px;
                border-left: 4px solid var(--vscode-textLink-foreground);
            }
            .stage h3 {
                margin: 0 0 5px 0;
                font-size: 16px;
            }
            .stage p {
                margin: 0;
                opacity: 0.8;
                font-size: 14px;
            }
        </style>
    </head>
    <body>
        <div class="dashboard-header">
            <h1>🔬 Overlay Science Lab Dashboard</h1>
            <p>Real-time monitoring of your scientific research pipeline</p>
        </div>
        
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">0</div>
                <div class="stat-label">Active Studies</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">0</div>
                <div class="stat-label">Deliverables Ready</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">6</div>
                <div class="stat-label">Pipeline Agents</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">7</div>
                <div class="stat-label">Scientific APIs</div>
            </div>
        </div>

        <div class="pipeline-stages">
            <h2>Pipeline Stages</h2>
            <div class="stage">
                <h3>📥 Stage 1: Source Ingestion</h3>
                <p>Parse PDFs, CSVs, and research materials</p>
            </div>
            <div class="stage">
                <h3>⚛️ Stage 2: Quantum Chaos Analysis</h3>
                <p>Run GLCCE quantum simulations</p>
            </div>
            <div class="stage">
                <h3>🧬 Stage 3: Biotech Simulation</h3>
                <p>Digital Lab protein and pathway analysis</p>
            </div>
            <div class="stage">
                <h3>📊 Stage 4: Insight Synthesis</h3>
                <p>DataLite analytics and visualization</p>
            </div>
            <div class="stage">
                <h3>📝 Stage 5: Paper & Book Writing</h3>
                <p>Generate research papers and book chapters</p>
            </div>
            <div class="stage">
                <h3>🎨 Stage 6: Metaphor Enhancement</h3>
                <p>Comic Engine visual metaphors</p>
            </div>
        </div>
    </body>
    </html>`;
}

/**
 * Create a new study
 */
async function createNewStudy() {
    const studyName = await vscode.window.showInputBox({
        prompt: 'Enter study name',
        placeHolder: 'e.g., Protein Folding Analysis 2026'
    });

    if (studyName) {
        const studyId = studyName.toLowerCase().replace(/\s+/g, '-');
        const config = vscode.workspace.getConfiguration('overlay-science-lab');
        const workspacePath = config.get('workspacePath', './science_workspace');
        const workspaceFolders = vscode.workspace.workspaceFolders;
        
        if (workspaceFolders) {
            const rootPath = workspaceFolders[0].uri.fsPath;
            const studyPath = path.join(rootPath, workspacePath, 'input', studyId);
            
            await fs.mkdir(studyPath, { recursive: true });
            vscode.window.showInformationMessage(`Study "${studyName}" created. Add source files to: ${studyPath}`);
            
            // Open the folder in explorer
            const uri = vscode.Uri.file(studyPath);
            await vscode.commands.executeCommand('revealFileInOS', uri);
        }
    }
}

/**
 * Upload sources to a study
 */
async function uploadSources() {
    const files = await vscode.window.showOpenDialog({
        canSelectMany: true,
        canSelectFiles: true,
        canSelectFolders: false,
        filters: {
            'Research Files': ['pdf', 'csv', 'json', 'txt', 'md', 'png', 'jpg']
        },
        title: 'Select Research Source Files'
    });

    if (files && files.length > 0) {
        vscode.window.showInformationMessage(`${files.length} source file(s) selected. Start pipeline?`, 'Yes', 'No')
            .then(answer => {
                if (answer === 'Yes') {
                    const studyId = `study_${Date.now()}`;
                    startPipeline(studyId);
                }
            });
    }
}

/**
 * View deliverables
 */
async function viewDeliverables() {
    const config = vscode.workspace.getConfiguration('overlay-science-lab');
    const workspacePath = config.get('workspacePath', './science_workspace');
    const workspaceFolders = vscode.workspace.workspaceFolders;
    
    if (workspaceFolders) {
        const rootPath = workspaceFolders[0].uri.fsPath;
        const deliverablesPath = path.join(rootPath, workspacePath, 'deliverables');
        const uri = vscode.Uri.file(deliverablesPath);
        await vscode.commands.executeCommand('revealFileInOS', uri);
    }
}

/**
 * Run quantum simulation
 */
async function runQuantumSimulation() {
    const chaosCoeff = await vscode.window.showInputBox({
        prompt: 'Enter chaos coefficient (0.0 - 1.0)',
        placeHolder: '0.7',
        value: '0.7'
    });

    if (chaosCoeff) {
        vscode.window.showInformationMessage(`Running GLCCE with chaos coefficient: ${chaosCoeff}`);
        statusBarItem.text = '$(sync~spin) Running Quantum Simulation...';
        
        setTimeout(() => {
            statusBarItem.text = '$(beaker) Science Lab: Ready';
            vscode.window.showInformationMessage('Quantum simulation complete!');
        }, 5000);
    }
}

/**
 * Run biotech simulation
 */
async function runBiotechSimulation() {
    const simType = await vscode.window.showQuickPick([
        'Protein Binding Study',
        'Pathway Analysis',
        'Molecular Dynamics',
        'Drug Screening'
    ], {
        placeHolder: 'Select simulation type'
    });

    if (simType) {
        vscode.window.showInformationMessage(`Running Digital Lab: ${simType}`);
        statusBarItem.text = '$(sync~spin) Running Biotech Simulation...';
        
        setTimeout(() => {
            statusBarItem.text = '$(beaker) Science Lab: Ready';
            vscode.window.showInformationMessage('Biotech simulation complete!');
        }, 5000);
    }
}

/**
 * Open science notebook
 */
async function openScienceNotebook() {
    const workspaceFolders = vscode.workspace.workspaceFolders;
    if (!workspaceFolders) {
        vscode.window.showErrorMessage('No workspace folder open');
        return;
    }

    const rootPath = workspaceFolders[0].uri.fsPath;
    const notebookPath = path.join(rootPath, 'science-notebook.science-nb');
    
    // Create a new notebook file if it doesn't exist
    try {
        await fs.access(notebookPath);
    } catch {
        await fs.writeFile(notebookPath, JSON.stringify({
            cells: [],
            metadata: {
                kernelspec: {
                    display_name: 'Python 3',
                    language: 'python',
                    name: 'python3'
                }
            }
        }, null, 2));
    }

    const uri = vscode.Uri.file(notebookPath);
    await vscode.commands.executeCommand('vscode.open', uri);
}

/**
 * Insert API snippet
 */
async function insertApiSnippet() {
    const api = await vscode.window.showQuickPick([
        'PubChemPy - Compound Search',
        'Biopython - PubMed Query',
        'BioThings - Gene Query',
        'PyBioPortal - Cancer Data',
        'ChEMBL - Drug Data',
        'Scikit-Bio - Sequence Analysis',
        'BioMart - Gene Attributes'
    ], {
        placeHolder: 'Select scientific API snippet'
    });

    if (api) {
        const editor = vscode.window.activeTextEditor;
        if (editor) {
            const snippetMap = {
                'PubChemPy - Compound Search': 'import pubchempy as pcp\ncompounds = pcp.get_compounds("aspirin", "name")\nfor compound in compounds:\n    print(compound.molecular_formula)',
                'Biopython - PubMed Query': 'from Bio import Entrez\nEntrez.email = "your@email.com"\nhandle = Entrez.esearch(db="pubmed", term="cancer research")\nrecord = Entrez.read(handle)\nprint(record["IdList"])',
                'BioThings - Gene Query': 'from biothings_client import get_client\nmg = get_client("gene")\nresults = mg.querymany(["TP53", "BRCA1"], scopes="symbol")\nfor result in results:\n    print(result)',
                'PyBioPortal - Cancer Data': 'from pybioportal import DataHub\ndh = DataHub()\nstudies = dh.get_all_studies()\nfor study in studies[:5]:\n    print(study)',
                'ChEMBL - Drug Data': 'from chembl_downloader import download\ndownload("molecule", "./chembl_data")',
                'Scikit-Bio - Sequence Analysis': 'from skbio import DNA\nseq = DNA("ACGTACGTACGT")\nprint(seq.complement())\nprint(seq.reverse_complement())',
                'BioMart - Gene Attributes': 'from biomart import Ensembl37\nserver = Ensembl37()\nattributes = ["ensembl_gene_id", "external_gene_name"]\nresults = server.attributes(attributes=attributes)'
            };

            const snippet = snippetMap[api];
            if (snippet) {
                editor.edit(editBuilder => {
                    editBuilder.insert(editor.selection.active, snippet);
                });
            }
        }
    }
}

/**
 * Stop active pipeline
 */
async function stopActivePipeline() {
    if (activePipelines.size === 0) {
        vscode.window.showInformationMessage('No active pipelines to stop');
        return;
    }

    const studyIds = Array.from(activePipelines.keys());
    const studyId = await vscode.window.showQuickPick(studyIds, {
        placeHolder: 'Select pipeline to stop'
    });

    if (studyId) {
        activePipelines.delete(studyId);
        vscode.window.showInformationMessage(`Pipeline ${studyId} stopped`);
        statusBarItem.text = '$(beaker) Science Lab: Ready';
    }
}

/**
 * Tree data providers
 */
class StudyTreeProvider {
    constructor() {
        this._onDidChangeTreeData = new vscode.EventEmitter();
        this.onDidChangeTreeData = this._onDidChangeTreeData.event;
    }

    refresh() {
        this._onDidChangeTreeData.fire();
    }

    getTreeItem(element) {
        return element;
    }

    async getChildren() {
        // Return active studies
        return [
            new vscode.TreeItem('No active studies', vscode.TreeItemCollapsibleState.None)
        ];
    }
}

class DeliverableTreeProvider {
    constructor() {
        this._onDidChangeTreeData = new vscode.EventEmitter();
        this.onDidChangeTreeData = this._onDidChangeTreeData.event;
    }

    refresh() {
        this._onDidChangeTreeData.fire();
    }

    getTreeItem(element) {
        return element;
    }

    async getChildren() {
        return [
            new vscode.TreeItem('No deliverables yet', vscode.TreeItemCollapsibleState.None)
        ];
    }
}

class ToolsTreeProvider {
    getTreeItem(element) {
        return element;
    }

    async getChildren() {
        const tools = [
            { name: 'GLCCE (Quantum Chaos)', icon: '⚛️' },
            { name: 'Digital Lab (Biotech)', icon: '🧬' },
            { name: 'Book Lab (Writing)', icon: '📝' },
            { name: 'Comic Engine (Metaphors)', icon: '🎨' },
            { name: 'DataLite (Analytics)', icon: '📊' },
            { name: 'PubChemPy API', icon: '🧪' },
            { name: 'Biopython API', icon: '🔬' }
        ];

        return tools.map(tool => {
            const item = new vscode.TreeItem(
                `${tool.icon} ${tool.name}`,
                vscode.TreeItemCollapsibleState.None
            );
            return item;
        });
    }
}

class AgentTreeProvider {
    getTreeItem(element) {
        return element;
    }

    async getChildren() {
        const agents = [
            { name: 'Source Ingester', status: 'ready' },
            { name: 'Quantum Runner', status: 'ready' },
            { name: 'Biotech Simulator', status: 'ready' },
            { name: 'Analytics Compiler', status: 'ready' },
            { name: 'Writer Agent', status: 'ready' },
            { name: 'Metaphor Agent', status: 'ready' }
        ];

        return agents.map(agent => {
            const item = new vscode.TreeItem(
                `${agent.status === 'ready' ? '✓' : '○'} ${agent.name}`,
                vscode.TreeItemCollapsibleState.None
            );
            item.description = agent.status;
            return item;
        });
    }
}

function deactivate() {
    if (statusBarItem) {
        statusBarItem.dispose();
    }
}

module.exports = {
    activate,
    deactivate
};
