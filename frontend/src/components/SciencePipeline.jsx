import React, { useState, useEffect, useRef } from 'react'

const STAGE_ICONS = {
  ingest: '📥',
  quantum: '⚛️',
  biotech: '🧬',
  analytics: '📊',
  writing: '📝',
  metaphor: '🎨',
}

const STATUS_COLORS = {
  idle: 'text-slate-400',
  queued: 'text-yellow-400',
  running: 'text-blue-400',
  success: 'text-emerald-400',
  error: 'text-red-400',
}

export default function SciencePipeline({ apiBase, agents }) {
  const [studyId, setStudyId] = useState('protein-folding-v1')
  const [studyType, setStudyType] = useState('quantum_biotech')
  const [execution, setExecution] = useState(null)
  const [stages, setStages] = useState([])
  const [loading, setLoading] = useState(false)
  const [outputs, setOutputs] = useState([])
  const wsRef = useRef(null)

  useEffect(() => {
    fetch(`${apiBase}/pipeline/stages`)
      .then(r => r.json())
      .then(d => setStages(d.stages || []))
      .catch(console.error)
  }, [apiBase])

  // Close the WebSocket when the component unmounts to prevent stale connections
  useEffect(() => {
    return () => {
      if (wsRef.current) {
        wsRef.current.close()
        wsRef.current = null
      }
    }
  }, [])

  const startPipeline = async () => {
    setLoading(true)
    try {
      const res = await fetch(`${apiBase}/execute`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ study_id: studyId, study_type: studyType }),
      })
      const data = await res.json()
      setExecution(data)
      connectWs(data.execution_id)
    } catch (e) {
      console.error(e)
    } finally {
      setLoading(false)
    }
  }

  const connectWs = (execId) => {
    if (wsRef.current) wsRef.current.close()
    const wsUrl = apiBase.replace(/^https?/, (p) => (p === 'https' ? 'wss' : 'ws')).replace('/api/v1', '')
    const ws = new WebSocket(`${wsUrl}/ws/executions/${execId}`)
    ws.onmessage = (evt) => {
      const data = JSON.parse(evt.data)
      setExecution(data)
      if (data.outputs) setOutputs(data.outputs)
      if (data.status === 'success' || data.status === 'error') ws.close()
    }
    wsRef.current = ws
  }

  const stageProgress = (stageId) => execution?.stages?.[stageId] || null

  return (
    <div className="space-y-6">
      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        {/* Config panel */}
        <div className="bg-slate-800 rounded-xl p-6 border border-slate-700">
          <h2 className="text-lg font-semibold mb-4 text-indigo-300">⚙️ Study Configuration</h2>
          <div className="space-y-4">
            <div>
              <label className="block text-sm text-slate-400 mb-1">Study ID</label>
              <input
                className="w-full bg-slate-700 border border-slate-600 rounded-lg px-3 py-2 text-sm text-slate-100 focus:outline-none focus:border-indigo-500"
                value={studyId}
                onChange={e => setStudyId(e.target.value)}
              />
            </div>
            <div>
              <label className="block text-sm text-slate-400 mb-1">Study Type</label>
              <select
                className="w-full bg-slate-700 border border-slate-600 rounded-lg px-3 py-2 text-sm text-slate-100 focus:outline-none focus:border-indigo-500"
                value={studyType}
                onChange={e => setStudyType(e.target.value)}
              >
                <option value="quantum_biotech">Quantum + Biotech</option>
                <option value="molecular_dynamics">Molecular Dynamics</option>
                <option value="protein_folding">Protein Folding</option>
                <option value="drug_discovery">Drug Discovery</option>
              </select>
            </div>
            <button
              onClick={startPipeline}
              disabled={loading || execution?.status === 'running'}
              className="w-full bg-indigo-600 hover:bg-indigo-500 disabled:bg-slate-600 disabled:cursor-not-allowed text-white font-medium py-2.5 rounded-lg transition-colors flex items-center justify-center gap-2"
            >
              {loading ? '⏳ Starting…' : execution?.status === 'running' ? '⚡ Running…' : '▶ Run Pipeline'}
            </button>
          </div>

          {/* Active agents */}
          <div className="mt-6">
            <h3 className="text-sm font-semibold text-slate-400 mb-3">AGENTS</h3>
            <div className="space-y-2">
              {agents.slice(0, 5).map(a => (
                <div key={a.id} className="flex items-center justify-between text-xs">
                  <span className="text-slate-300">{a.name}</span>
                  <span className={a.status === 'busy' ? 'text-blue-400 animate-pulse' : 'text-emerald-400'}>
                    {a.status === 'busy' ? '⚡ busy' : '✓ idle'}
                  </span>
                </div>
              ))}
            </div>
          </div>
        </div>

        {/* Pipeline stages */}
        <div className="lg:col-span-2 bg-slate-800 rounded-xl p-6 border border-slate-700">
          <h2 className="text-lg font-semibold mb-4 text-indigo-300">🔄 Pipeline Stages</h2>

          {execution && (
            <div className="mb-4">
              <div className="flex items-center justify-between text-sm mb-1">
                <span className={STATUS_COLORS[execution.status] || 'text-slate-400'}>
                  {execution.message}
                </span>
                <span className="text-slate-400">{execution.progress}%</span>
              </div>
              <div className="w-full bg-slate-700 rounded-full h-2">
                <div
                  className="bg-indigo-500 h-2 rounded-full transition-all duration-500"
                  style={{ width: `${execution.progress}%` }}
                />
              </div>
            </div>
          )}

          <div className="space-y-3">
            {stages.map(stage => {
              const sp = stageProgress(stage.id)
              return (
                <div key={stage.id} className="flex items-center gap-4 bg-slate-700/50 rounded-lg px-4 py-3">
                  <span className="text-2xl">{STAGE_ICONS[stage.id] || '🔬'}</span>
                  <div className="flex-1 min-w-0">
                    <div className="flex items-center justify-between">
                      <span className="text-sm font-medium text-slate-200">{stage.name}</span>
                      <span className="text-xs text-slate-500">Order {stage.order}</span>
                    </div>
                    {sp && (
                      <div className="mt-1 w-full bg-slate-600 rounded-full h-1.5">
                        <div
                          className={`h-1.5 rounded-full transition-all duration-300 ${sp.status === 'complete' ? 'bg-emerald-500' : 'bg-blue-500'}`}
                          style={{ width: `${sp.progress}%` }}
                        />
                      </div>
                    )}
                  </div>
                  <span className="text-xs shrink-0">
                    {!sp && <span className="text-slate-500">pending</span>}
                    {sp?.status === 'running' && <span className="text-blue-400 animate-pulse">running</span>}
                    {sp?.status === 'complete' && <span className="text-emerald-400">✓ done</span>}
                  </span>
                </div>
              )
            })}
          </div>
        </div>
      </div>

      {/* Outputs */}
      {outputs.length > 0 && (
        <div className="bg-slate-800 rounded-xl p-6 border border-slate-700">
          <h2 className="text-lg font-semibold mb-4 text-emerald-300">📦 Generated Outputs</h2>
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4">
            {outputs.map(out => (
              <div key={out.id} className="bg-slate-700 rounded-lg p-4 border border-slate-600">
                <div className="text-2xl mb-2">
                  {out.type === 'paper' ? '📄' : out.type === 'book_chapter' ? '📚' : out.type === 'analytics' ? '📊' : '📁'}
                </div>
                <p className="text-sm font-medium text-slate-200 leading-tight">{out.title}</p>
                <div className="mt-2 flex flex-wrap gap-1">
                  {out.pages && <span className="text-xs bg-slate-600 px-2 py-0.5 rounded-full text-slate-300">{out.pages}p</span>}
                  {out.citations && <span className="text-xs bg-slate-600 px-2 py-0.5 rounded-full text-slate-300">{out.citations} refs</span>}
                  {out.charts && <span className="text-xs bg-slate-600 px-2 py-0.5 rounded-full text-slate-300">{out.charts} charts</span>}
                </div>
                <span className="mt-2 inline-block text-xs text-emerald-400">✓ {out.status}</span>
              </div>
            ))}
          </div>
        </div>
      )}
    </div>
  )
}
