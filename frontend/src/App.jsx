import React, { useState, useEffect } from 'react'
import SciencePipeline from './components/SciencePipeline'
import CrossDisciplinaryAgents from './components/CrossDisciplinaryAgents'
import CISAssistant from './components/CISAssistant'
import CirculatoryInformatics from './components/CirculatoryInformatics'

const API_BASE = 'http://localhost:8000/api/v1'

const TABS = [
  { id: 'pipeline', label: '🔬 Science Pipeline' },
  { id: 'agents', label: '🤖 Agents' },
  { id: 'cis', label: '🩺 CIS Assistant' },
  { id: 'circulatory', label: '🫀 Circulatory Informatics' },
]

export default function App() {
  const [activeTab, setActiveTab] = useState('pipeline')
  const [agents, setAgents] = useState([])
  const [disciplines, setDisciplines] = useState([])
  const [health, setHealth] = useState(null)

  useEffect(() => {
    Promise.all([
      fetch(`${API_BASE}/agents`).then(r => r.json()).then(d => setAgents(d.agents || [])),
      fetch(`${API_BASE}/disciplines`).then(r => r.json()).then(d => setDisciplines(d.disciplines || [])),
      fetch(`${API_BASE}/health`).then(r => r.json()).then(setHealth),
    ]).catch(console.error)
  }, [])

  return (
    <div className="min-h-screen bg-slate-900 text-slate-100">
      {/* Header */}
      <header className="bg-slate-800 border-b border-slate-700 px-6 py-4">
        <div className="flex items-center justify-between max-w-7xl mx-auto">
          <div>
            <h1 className="text-2xl font-bold text-indigo-400">🧬 Overlay Science Team</h1>
            <p className="text-slate-400 text-sm mt-0.5">Agentic Cross-Disciplinary Research Factory</p>
          </div>
          <div className="flex items-center gap-3">
            {health && (
              <span className="flex items-center gap-1.5 text-sm text-emerald-400">
                <span className="w-2 h-2 rounded-full bg-emerald-400 animate-pulse" />
                {health.status} · {health.agents_active} agents active
              </span>
            )}
          </div>
        </div>
      </header>

      {/* Navigation */}
      <nav className="bg-slate-800 border-b border-slate-700 px-6">
        <div className="flex max-w-7xl mx-auto">
          {TABS.map(tab => (
            <button
              key={tab.id}
              onClick={() => setActiveTab(tab.id)}
              className={`px-5 py-3 text-sm font-medium transition-colors border-b-2 ${
                activeTab === tab.id
                  ? 'border-indigo-400 text-indigo-400'
                  : 'border-transparent text-slate-400 hover:text-slate-200'
              }`}
            >
              {tab.label}
            </button>
          ))}
        </div>
      </nav>

      {/* Content */}
      <main className="max-w-7xl mx-auto px-6 py-8">
        {activeTab === 'pipeline' && <SciencePipeline apiBase={API_BASE} agents={agents} />}
        {activeTab === 'agents' && <CrossDisciplinaryAgents agents={agents} disciplines={disciplines} />}
        {activeTab === 'cis' && <CISAssistant apiBase={API_BASE} />}
        {activeTab === 'circulatory' && <CirculatoryInformatics agents={agents} health={health} />}
      </main>
    </div>
  )
}
