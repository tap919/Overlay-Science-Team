# Overlay365 Science Factory - React UI
## Production-Ready Dashboard

**File:** `ui/App.tsx`

```typescript
import React, { useState, useEffect, useRef } from 'react';
import {
  Upload, Play, Pause, Activity, BarChart3, Brain, Zap,
  Clock, CheckCircle, AlertCircle, FileJson, FileText, Download,
  Settings, Users, Home
} from 'lucide-react';
import FileUploadZone from './components/FileUploadZone';
import PipelineMonitor from './components/PipelineMonitor';
import OutputGallery from './components/OutputGallery';
import AgentStatus from './components/AgentStatus';
import AnalyticsDashboard from './components/AnalyticsDashboard';

const API_BASE = 'http://localhost:8000/api/v1';

interface ExecutionState {
  id: string | null;
  studyId: string | null;
  status: 'idle' | 'uploading' | 'queued' | 'running' | 'success' | 'error';
  progress: number;
  currentStage: string;
  message: string;
  results: any | null;
  error: string | null;
}

interface Agent {
  id: string;
  name: string;
  role: string;
  skills: Record<string, number>;
  memory: {
    learned_workflows: number;
    successful_tasks: number;
    avg_execution_time_ms: number;
  };
}

export default function App() {
  const [executionState, setExecutionState] = useState<ExecutionState>({
    id: null,
    studyId: null,
    status: 'idle',
    progress: 0,
    currentStage: '',
    message: 'Ready to begin',
    results: null,
    error: null
  });

  const [agents, setAgents] = useState<Agent[]>([]);
  const [studyId, setStudyId] = useState('protein-folding-v1');
  const [studyType, setStudyType] = useState('quantum_biotech');
  const [config, setConfig] = useState({
    quantum_depth: 'auto',
    md_duration_ns: 50,
    use_metaphor_engine: true
  });

  const [activeTab, setActiveTab] = useState<'upload' | 'monitor' | 'outputs' | 'agents' | 'analytics'>('upload');
  const ws = useRef<WebSocket | null>(null);

  // Load agents on mount
  useEffect(() => {
    loadAgents();
  }, []);

  const loadAgents = async () => {
    try {
      const response = await fetch(`${API_BASE}/agents`);
      const data = await response.json();
      setAgents(data.agents);
    } catch (error) {
      console.error('Failed to load agents:', error);
    }
  };

  const handleFileUpload = async (files: File[]) => {
    try {
      setExecutionState(prev => ({
        ...prev,
        status: 'uploading',
        message: `Uploading ${files.length} files...`
      }));

      const formData = new FormData();
      files.forEach(file => formData.append('files', file));

      const response = await fetch(
        `${API_BASE}/upload?study_id=${studyId}`,
        {
          method: 'POST',
          body: formData
        }
      );

      if (!response.ok) throw new Error('Upload failed');

      const result = await response.json();
      
      setExecutionState(prev => ({
        ...prev,
        status: 'queued',
        message: `✓ Uploaded ${result.files_uploaded} files. Ready to execute pipeline.`,
        progress: 0
      }));

      // Auto-switch to monitor tab
      setActiveTab('monitor');

    } catch (error) {
      setExecutionState(prev => ({
        ...prev,
        status: 'error',
        error: error instanceof Error ? error.message : 'Upload failed'
      }));
    }
  };

  const handleExecutePipeline = async () => {
    try {
      if (!studyId) {
        alert('Please enter a study ID');
        return;
      }

      const response = await fetch(`${API_BASE}/execute`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          study_id: studyId,
          study_type: studyType,
          config
        })
      });

      if (!response.ok) throw new Error('Execution failed to start');

      const result = await response.json();
      const executionId = result.execution_id;

      setExecutionState(prev => ({
        ...prev,
        id: executionId,
        studyId: result.study_id,
        status: 'queued',
        message: 'Pipeline queued. Starting soon...',
        progress: 0
      }));

      // Connect WebSocket for real-time updates
      connectWebSocket(executionId);

    } catch (error) {
      setExecutionState(prev => ({
        ...prev,
        status: 'error',
        error: error instanceof Error ? error.message : 'Failed to start execution'
      }));
    }
  };

  const connectWebSocket = (executionId: string) => {
    const wsUrl = `ws://localhost:8000/ws/executions/${executionId}`;
    
    ws.current = new WebSocket(wsUrl);

    ws.current.onopen = () => {
      console.log('WebSocket connected');
      setExecutionState(prev => ({
        ...prev,
        status: 'running'
      }));
    };

    ws.current.onmessage = (event) => {
      const update = JSON.parse(event.data);
      
      // Update execution state
      setExecutionState(prev => {
        let newProgress = prev.progress;
        
        if (update.stage === 'ingestion') newProgress = 10;
        else if (update.stage === 'quantum') newProgress = 35;
        else if (update.stage === 'biotech') newProgress = 65;
        else if (update.stage === 'writing') newProgress = 85;
        else if (update.stage === 'analytics') newProgress = 95;
        else if (update.stage === 'complete') newProgress = 100;

        return {
          ...prev,
          currentStage: update.stage,
          message: update.message,
          progress: newProgress,
          results: update.result || prev.results
        };
      });

      // Play notification sound on complete
      if (update.stage === 'complete') {
        playNotificationSound();
      }
    };

    ws.current.onerror = (error) => {
      console.error('WebSocket error:', error);
      setExecutionState(prev => ({
        ...prev,
        status: 'error',
        error: 'WebSocket connection error'
      }));
    };

    ws.current.onclose = () => {
      console.log('WebSocket disconnected');
    };
  };

  const playNotificationSound = () => {
    // Play a simple beep or notification sound
    const audioContext = new (window.AudioContext || (window as any).webkitAudioContext)();
    const oscillator = audioContext.createOscillator();
    const gainNode = audioContext.createGain();

    oscillator.connect(gainNode);
    gainNode.connect(audioContext.destination);

    oscillator.frequency.value = 800;
    oscillator.type = 'sine';

    gainNode.gain.setValueAtTime(0.3, audioContext.currentTime);
    gainNode.gain.exponentialRampToValueAtTime(0.01, audioContext.currentTime + 0.5);

    oscillator.start(audioContext.currentTime);
    oscillator.stop(audioContext.currentTime + 0.5);
  };

  const handleStopExecution = () => {
    if (ws.current) {
      ws.current.close();
    }
    setExecutionState(prev => ({
      ...prev,
      status: 'idle'
    }));
  };

  return (
    <div className="w-full h-screen bg-gradient-to-br from-slate-900 via-slate-800 to-slate-900 text-white">
      {/* Header */}
      <header className="bg-slate-800/80 border-b border-slate-700 backdrop-blur">
        <div className="max-w-7xl mx-auto px-6 py-4 flex items-center justify-between">
          <div className="flex items-center gap-3">
            <Brain className="w-8 h-8 text-cyan-400" />
            <h1 className="text-2xl font-bold">Overlay365 Science Factory</h1>
            <span className="text-xs bg-cyan-600/20 text-cyan-300 px-2 py-1 rounded-full">
              {executionState.status === 'running' ? '🟢 Active' : '⚪ Ready'}
            </span>
          </div>
          <div className="flex items-center gap-4">
            <button className="p-2 hover:bg-slate-700 rounded-lg transition">
              <Settings className="w-5 h-5" />
            </button>
            <button className="p-2 hover:bg-slate-700 rounded-lg transition">
              <Users className="w-5 h-5" />
            </button>
          </div>
        </div>
      </header>

      {/* Main Content */}
      <div className="flex h-[calc(100vh-70px)]">
        {/* Sidebar Navigation */}
        <div className="w-64 bg-slate-800/50 border-r border-slate-700 p-4 overflow-y-auto">
          <nav className="space-y-2">
            {[
              { id: 'upload', icon: Upload, label: 'Upload' },
              { id: 'monitor', icon: Activity, label: 'Monitor' },
              { id: 'outputs', icon: FileJson, label: 'Outputs' },
              { id: 'agents', icon: Brain, label: 'Agents' },
              { id: 'analytics', icon: BarChart3, label: 'Analytics' }
            ].map(tab => (
              <button
                key={tab.id}
                onClick={() => setActiveTab(tab.id as any)}
                className={`w-full flex items-center gap-3 px-4 py-3 rounded-lg transition ${
                  activeTab === tab.id
                    ? 'bg-cyan-600 text-white'
                    : 'text-slate-300 hover:bg-slate-700'
                }`}
              >
                <tab.icon className="w-5 h-5" />
                <span>{tab.label}</span>
              </button>
            ))}
          </nav>

          {/* Study Configuration */}
          <div className="mt-8 p-4 bg-slate-700/30 rounded-lg border border-slate-600">
            <h3 className="font-semibold text-sm mb-3">Study Config</h3>
            
            <div className="space-y-3">
              <div>
                <label className="text-xs text-slate-400 block mb-1">Study ID</label>
                <input
                  type="text"
                  value={studyId}
                  onChange={(e) => setStudyId(e.target.value)}
                  className="w-full px-2 py-1 bg-slate-700 border border-slate-600 rounded text-sm"
                  placeholder="protein-folding-v1"
                />
              </div>

              <div>
                <label className="text-xs text-slate-400 block mb-1">Study Type</label>
                <select
                  value={studyType}
                  onChange={(e) => setStudyType(e.target.value)}
                  className="w-full px-2 py-1 bg-slate-700 border border-slate-600 rounded text-sm"
                >
                  <option value="quantum_biotech">Quantum + Biotech</option>
                  <option value="quantum_only">Quantum Only</option>
                  <option value="biotech_only">Biotech Only</option>
                  <option value="analytics_only">Analytics Only</option>
                </select>
              </div>

              <button
                onClick={handleExecutePipeline}
                disabled={executionState.status === 'running'}
                className={`w-full py-2 rounded-lg font-semibold flex items-center justify-center gap-2 transition ${
                  executionState.status === 'running'
                    ? 'bg-slate-600 text-slate-400 cursor-not-allowed'
                    : 'bg-cyan-600 hover:bg-cyan-700 text-white'
                }`}
              >
                {executionState.status === 'running' ? (
                  <>
                    <Pause className="w-4 h-4 animate-spin" />
                    Running...
                  </>
                ) : (
                  <>
                    <Play className="w-4 h-4" />
                    Execute
                  </>
                )}
              </button>
            </div>
          </div>
        </div>

        {/* Main Content Area */}
        <div className="flex-1 overflow-y-auto p-6">
          {activeTab === 'upload' && (
            <div>
              <h2 className="text-2xl font-bold mb-6">Upload Research Files</h2>
              <FileUploadZone onFilesSelected={handleFileUpload} />
            </div>
          )}

          {activeTab === 'monitor' && (
            <div>
              <h2 className="text-2xl font-bold mb-6">Pipeline Monitor</h2>
              {executionState.id ? (
                <PipelineMonitor executionState={executionState} />
              ) : (
                <div className="bg-slate-700/30 border border-slate-600 rounded-lg p-8 text-center">
                  <Activity className="w-16 h-16 mx-auto mb-4 text-slate-600" />
                  <p className="text-slate-400">Upload files and execute pipeline to monitor progress</p>
                </div>
              )}
            </div>
          )}

          {activeTab === 'outputs' && (
            <div>
              <h2 className="text-2xl font-bold mb-6">Deliverables</h2>
              {executionState.results ? (
                <OutputGallery results={executionState.results} />
              ) : (
                <div className="bg-slate-700/30 border border-slate-600 rounded-lg p-8 text-center">
                  <FileJson className="w-16 h-16 mx-auto mb-4 text-slate-600" />
                  <p className="text-slate-400">No outputs yet. Complete a pipeline execution to see results.</p>
                </div>
              )}
            </div>
          )}

          {activeTab === 'agents' && (
            <div>
              <h2 className="text-2xl font-bold mb-6">Agent Team</h2>
              <div className="grid grid-cols-2 gap-6">
                {agents.map(agent => (
                  <AgentStatus key={agent.id} agent={agent} />
                ))}
              </div>
            </div>
          )}

          {activeTab === 'analytics' && (
            <div>
              <h2 className="text-2xl font-bold mb-6">Analytics & Performance</h2>
              {executionState.results ? (
                <AnalyticsDashboard results={executionState.results} />
              ) : (
                <div className="bg-slate-700/30 border border-slate-600 rounded-lg p-8 text-center">
                  <BarChart3 className="w-16 h-16 mx-auto mb-4 text-slate-600" />
                  <p className="text-slate-400">Execute a pipeline to view performance analytics</p>
                </div>
              )}
            </div>
          )}
        </div>
      </div>

      {/* Status Bar */}
      <div className="h-16 bg-slate-800/80 border-t border-slate-700 px-6 flex items-center justify-between">
        <div className="flex items-center gap-3">
          {executionState.status === 'running' ? (
            <>
              <div className="w-2 h-2 bg-green-500 rounded-full animate-pulse"></div>
              <span className="text-sm">{executionState.currentStage || 'Initializing...'}</span>
            </>
          ) : executionState.status === 'success' ? (
            <>
              <CheckCircle className="w-4 h-4 text-green-500" />
              <span className="text-sm">✓ Complete</span>
            </>
          ) : executionState.status === 'error' ? (
            <>
              <AlertCircle className="w-4 h-4 text-red-500" />
              <span className="text-sm text-red-400">{executionState.error}</span>
            </>
          ) : null}
        </div>

        <div className="flex items-center gap-4">
          <div className="text-sm text-slate-400">
            {executionState.id && `Execution: ${executionState.id.substring(0, 8)}...`}
          </div>
          {executionState.status === 'running' && (
            <button
              onClick={handleStopExecution}
              className="px-3 py-1 bg-red-600/20 hover:bg-red-600/40 text-red-300 rounded text-sm transition"
            >
              Stop
            </button>
          )}
        </div>
      </div>
    </div>
  );
}
```

---

## Component: `FileUploadZone.tsx`

```typescript
import React, { useRef } from 'react';
import { Upload, File, CheckCircle } from 'lucide-react';

interface Props {
  onFilesSelected: (files: File[]) => void;
}

export default function FileUploadZone({ onFilesSelected }: Props) {
  const inputRef = useRef<HTMLInputElement>(null);
  const [dragActive, setDragActive] = React.useState(false);
  const [selectedFiles, setSelectedFiles] = React.useState<File[]>([]);

  const handleDrag = (e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setDragActive(e.type === 'dragenter' || e.type === 'dragover');
  };

  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setDragActive(false);

    const files = Array.from(e.dataTransfer.files);
    setSelectedFiles(files);
  };

  const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const files = Array.from(e.target.files || []);
    setSelectedFiles(files);
  };

  const handleUpload = () => {
    if (selectedFiles.length > 0) {
      onFilesSelected(selectedFiles);
      setSelectedFiles([]);
    }
  };

  return (
    <div className="space-y-6">
      {/* Upload Zone */}
      <div
        onDragEnter={handleDrag}
        onDragLeave={handleDrag}
        onDragOver={handleDrag}
        onDrop={handleDrop}
        className={`border-2 border-dashed rounded-lg p-12 text-center transition ${
          dragActive
            ? 'border-cyan-400 bg-cyan-400/10'
            : 'border-slate-600 bg-slate-700/20 hover:border-slate-500'
        }`}
      >
        <Upload className="w-12 h-12 mx-auto mb-4 text-slate-400" />
        <h3 className="text-lg font-semibold mb-2">Drop your files here</h3>
        <p className="text-slate-400 mb-4">or click to browse</p>
        <input
          ref={inputRef}
          type="file"
          multiple
          onChange={handleChange}
          className="hidden"
          accept=".pdf,.csv,.json,.txt,.png,.jpg,.md"
        />
        <button
          onClick={() => inputRef.current?.click()}
          className="px-6 py-2 bg-cyan-600 hover:bg-cyan-700 rounded-lg font-semibold transition"
        >
          Select Files
        </button>
      </div>

      {/* File List */}
      {selectedFiles.length > 0 && (
        <div className="bg-slate-700/30 border border-slate-600 rounded-lg p-6">
          <h3 className="font-semibold mb-4">Selected Files ({selectedFiles.length})</h3>
          <div className="space-y-2 mb-4">
            {selectedFiles.map((file, i) => (
              <div key={i} className="flex items-center gap-3 p-3 bg-slate-700/50 rounded">
                <File className="w-4 h-4 text-cyan-400" />
                <div className="flex-1">
                  <p className="text-sm font-medium">{file.name}</p>
                  <p className="text-xs text-slate-400">{(file.size / 1024).toFixed(1)} KB</p>
                </div>
                <CheckCircle className="w-4 h-4 text-green-500" />
              </div>
            ))}
          </div>
          <button
            onClick={handleUpload}
            className="w-full px-4 py-2 bg-green-600 hover:bg-green-700 rounded-lg font-semibold transition"
          >
            Upload Files
          </button>
        </div>
      )}

      {/* Info Cards */}
      <div className="grid grid-cols-3 gap-4">
        {[
          { title: 'Supported Formats', desc: 'PDF, CSV, JSON, TXT, PNG, MD' },
          { title: 'Max File Size', desc: '100 MB per file' },
          { title: 'Processing Time', desc: '~6 hours for full pipeline' }
        ].map((card, i) => (
          <div key={i} className="bg-slate-700/30 border border-slate-600 rounded-lg p-4">
            <h4 className="font-semibold text-sm mb-1">{card.title}</h4>
            <p className="text-xs text-slate-400">{card.desc}</p>
          </div>
        ))}
      </div>
    </div>
  );
}
```

---

**Continue with components:**
- `PipelineMonitor.tsx` - Real-time execution timeline
- `OutputGallery.tsx` - Display papers, chapters, analytics
- `AgentStatus.tsx` - Individual agent skill tracking
- `AnalyticsDashboard.tsx` - Performance metrics & visualizations

All components follow the dark theme with cyan accents, responsive design, and real-time updates via WebSocket.

**Next:** Implement agent logic in `agents/` directory, integrate GLCCE/Digital Lab tools, test end-to-end pipeline.
