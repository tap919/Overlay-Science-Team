const vscode = require('vscode');
const path = require('path');
const fs = require('fs').promises;
const fetch = require('node-fetch');

const API_BASE = 'http://localhost:8000/api/v1';

/**
 * End-to-end agentic automation orchestrator.
 * - Collects sources from workspace
 * - Registers study + uploads sources to backend
 * - Kicks off pipeline execution
 * - Streams status to status bar and refreshes tree views
 */
async function runAgenticWorkflow({ statusBarItem, studyTreeProvider, deliverableTreeProvider }) {
    const studyId = await vscode.window.showInputBox({
        prompt: 'Enter study ID for agentic automation',
        placeHolder: 'e.g., auto-study-protein-folding',
        value: `auto-study-${Date.now()}`
    });
    if (!studyId) {
        return;
    }

    // Locate sources under configured workspace
    const config = vscode.workspace.getConfiguration('overlay-science-lab');
    const workspacePath = config.get('workspacePath', './science_workspace');
    const workspaceFolders = vscode.workspace.workspaceFolders;
    if (!workspaceFolders || workspaceFolders.length === 0) {
        vscode.window.showErrorMessage('No workspace folder open. Open a folder to run automation.');
        return;
    }
    const rootPath = workspaceFolders[0].uri.fsPath;
    const sourceDir = path.join(rootPath, workspacePath, 'input');

    let files = [];
    try {
        const entries = await fs.readdir(sourceDir, { withFileTypes: true });
        for (const entry of entries) {
            if (entry.isFile()) {
                files.push(path.join(sourceDir, entry.name));
            }
        }
    } catch (err) {
        vscode.window.showErrorMessage(`Unable to read sources from ${sourceDir}`);
        return;
    }

    if (files.length === 0) {
        vscode.window.showWarningMessage(`No source files found in ${sourceDir}. Add files then retry.`);
        return;
    }

    statusBarItem.text = '$(sync~spin) Agentic workflow: uploading sources...';

    try {
        const formData = new (require('form-data'))();
        files.forEach(f => {
            formData.append('files', require('fs').createReadStream(f));
        });

        const uploadResp = await fetch(`${API_BASE}/upload?study_id=${encodeURIComponent(studyId)}`, {
            method: 'POST',
            body: formData
        });
        if (!uploadResp.ok) {
            throw new Error(`Upload failed: ${uploadResp.statusText}`);
        }

        statusBarItem.text = '$(sync~spin) Agentic workflow: starting pipeline...';
        const execResp = await fetch(`${API_BASE}/execute`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                study_id: studyId,
                study_type: 'quantum_biotech',
                config: {
                    quantum_depth: 'auto',
                    md_duration_ns: 50,
                    use_metaphor_engine: true,
                    batch_mode: true
                }
            })
        });
        if (!execResp.ok) {
            throw new Error(`Execution failed: ${execResp.statusText}`);
        }
        const exec = await execResp.json();
        const executionId = exec.execution_id;

        statusBarItem.text = '$(sync~spin) Agentic workflow: monitoring pipeline...';
        await monitorExecution(executionId, statusBarItem, deliverableTreeProvider);
        studyTreeProvider.refresh();
    } catch (err) {
        statusBarItem.text = '$(alert) Agentic workflow failed';
        vscode.window.showErrorMessage(err.message);
    }
}

async function monitorExecution(executionId, statusBarItem, deliverableTreeProvider) {
    const wsUrl = `ws://localhost:8000/ws/executions/${executionId}`;
    return new Promise((resolve) => {
        const ws = new (require('ws'))(wsUrl);

        ws.on('open', () => {
            statusBarItem.text = '$(sync~spin) Agentic workflow: running...';
        });

        ws.on('message', (data) => {
            try {
                const update = JSON.parse(data.toString());
                if (update.stage) {
                    statusBarItem.text = `$(sync~spin) ${update.stage}: ${update.message || ''}`;
                }
                if (update.stage === 'complete') {
                    statusBarItem.text = '$(check) Agentic workflow complete';
                    deliverableTreeProvider.refresh();
                    ws.close();
                    resolve();
                } else if (update.stage === 'error') {
                    statusBarItem.text = '$(alert) Agentic workflow error';
                    ws.close();
                    resolve();
                }
            } catch (e) {
                // ignore malformed updates
            }
        });

        ws.on('error', () => {
            statusBarItem.text = '$(alert) Agentic workflow connection error';
            resolve();
        });
    });
}

module.exports = {
    runAgenticWorkflow
};
