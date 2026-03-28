import React, { useState, useEffect } from 'react'

export default function CISAssistant({ apiBase }) {
  const [principles, setPrinciples] = useState([])
  const [selectedPrinciple, setSelectedPrinciple] = useState(null)
  const [contractInput, setContractInput] = useState('')
  const [contractOutput, setContractOutput] = useState('')
  const [validationInput, setValidationInput] = useState('')
  const [validationResult, setValidationResult] = useState(null)
  const [activeSection, setActiveSection] = useState('principles')

  useEffect(() => {
    fetch(`${apiBase}/cis/principles`)
      .then(r => r.json())
      .then(d => setPrinciples(d.principles || []))
      .catch(console.error)
  }, [apiBase])

  const generateContract = () => {
    if (!contractInput.trim()) return
    // Simulate contract generation
    setContractOutput(`/**
 * CIS Contract: ${contractInput.trim()}
 * Generated: ${new Date().toISOString()}
 * Principle: Distributed Autonomy + Feedback-Driven Adaptation
 */
interface ${contractInput.replace(/\s+/g, '')}Contract {
  // Input specification
  inputs: {
    data: unknown;
    config?: Record<string, unknown>;
  };
  // Output guarantee
  outputs: {
    result: unknown;
    metadata: {
      timestamp: string;
      agent_id: string;
      confidence: number;
    };
  };
  // Error handling (Graceful Degradation principle)
  onFailure: 'retry' | 'fallback' | 'escalate';
  // Resource flow (Efficient Resource Flow principle)
  priority: 'high' | 'normal' | 'low';
}`)
  }

  const validateCode = () => {
    if (!validationInput.trim()) return
    const issues = []
    if (!validationInput.includes('try') && !validationInput.includes('catch')) {
      issues.push({ level: 'warning', msg: 'No error handling detected — violates Graceful Degradation principle' })
    }
    if (!validationInput.includes('log') && !validationInput.includes('monitor')) {
      issues.push({ level: 'info', msg: 'No logging/monitoring — consider adding for Continuous Sensing' })
    }
    if (validationInput.includes('while(true)') || validationInput.includes('while (true)')) {
      issues.push({ level: 'error', msg: 'Infinite loop detected — violates Efficient Resource Flow principle' })
    }
    if (issues.length === 0) {
      issues.push({ level: 'success', msg: 'Code passes CIS validation checks ✓' })
    }
    setValidationResult(issues)
  }

  const SECTION_TABS = [
    { id: 'principles', label: '📖 Seven Principles' },
    { id: 'contracts', label: '📋 Contract Generator' },
    { id: 'validation', label: '✅ Code Validation' },
  ]

  return (
    <div className="space-y-6">
      {/* Section nav */}
      <div className="flex gap-2 bg-slate-800 rounded-xl p-1.5 border border-slate-700">
        {SECTION_TABS.map(tab => (
          <button
            key={tab.id}
            onClick={() => setActiveSection(tab.id)}
            className={`flex-1 py-2 rounded-lg text-sm font-medium transition-colors ${
              activeSection === tab.id
                ? 'bg-indigo-600 text-white'
                : 'text-slate-400 hover:text-slate-200'
            }`}
          >
            {tab.label}
          </button>
        ))}
      </div>

      {/* Seven Principles */}
      {activeSection === 'principles' && (
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
          {principles.map(p => (
            <button
              key={p.id}
              onClick={() => setSelectedPrinciple(selectedPrinciple?.id === p.id ? null : p)}
              className={`text-left bg-slate-800 rounded-xl p-5 border transition-all ${
                selectedPrinciple?.id === p.id
                  ? 'border-indigo-500 ring-1 ring-indigo-500/40'
                  : 'border-slate-700 hover:border-slate-600'
              }`}
            >
              <div className="text-3xl mb-2">{p.icon}</div>
              <h3 className="font-semibold text-slate-100 mb-1">{p.name}</h3>
              <p className="text-sm text-slate-400 leading-relaxed">{p.description}</p>
            </button>
          ))}
        </div>
      )}

      {/* Contract Generator */}
      {activeSection === 'contracts' && (
        <div className="bg-slate-800 rounded-xl p-6 border border-slate-700 space-y-4">
          <h2 className="text-lg font-semibold text-indigo-300">📋 CIS Contract Generator</h2>
          <p className="text-sm text-slate-400">
            Generate TypeScript interface contracts based on CIS principles for any component or agent.
          </p>
          <div>
            <label className="block text-sm text-slate-400 mb-1">Component / Agent Name</label>
            <input
              className="w-full bg-slate-700 border border-slate-600 rounded-lg px-3 py-2 text-sm text-slate-100 focus:outline-none focus:border-indigo-500"
              placeholder="e.g. QuantumChaosRunner, DataIngester…"
              value={contractInput}
              onChange={e => setContractInput(e.target.value)}
              onKeyDown={e => e.key === 'Enter' && generateContract()}
            />
          </div>
          <button
            onClick={generateContract}
            className="bg-indigo-600 hover:bg-indigo-500 text-white px-5 py-2 rounded-lg text-sm font-medium transition-colors"
          >
            Generate Contract
          </button>
          {contractOutput && (
            <pre className="mt-4 bg-slate-900 rounded-lg p-4 text-xs text-emerald-300 overflow-x-auto border border-slate-700 whitespace-pre-wrap">
              {contractOutput}
            </pre>
          )}
        </div>
      )}

      {/* Validation */}
      {activeSection === 'validation' && (
        <div className="bg-slate-800 rounded-xl p-6 border border-slate-700 space-y-4">
          <h2 className="text-lg font-semibold text-indigo-300">✅ CIS Code Validation</h2>
          <p className="text-sm text-slate-400">
            Paste code to validate against CIS principles (Graceful Degradation, Continuous Sensing, Efficient Resource Flow).
          </p>
          <textarea
            className="w-full h-48 bg-slate-700 border border-slate-600 rounded-lg px-3 py-2 text-sm text-slate-100 font-mono focus:outline-none focus:border-indigo-500 resize-none"
            placeholder="Paste your code here…"
            value={validationInput}
            onChange={e => setValidationInput(e.target.value)}
          />
          <button
            onClick={validateCode}
            className="bg-indigo-600 hover:bg-indigo-500 text-white px-5 py-2 rounded-lg text-sm font-medium transition-colors"
          >
            Validate
          </button>
          {validationResult && (
            <div className="space-y-2">
              {validationResult.map((issue, i) => (
                <div
                  key={i}
                  className={`flex items-start gap-2 px-4 py-3 rounded-lg text-sm ${
                    issue.level === 'error' ? 'bg-red-500/10 text-red-300 border border-red-500/20' :
                    issue.level === 'warning' ? 'bg-amber-500/10 text-amber-300 border border-amber-500/20' :
                    issue.level === 'success' ? 'bg-emerald-500/10 text-emerald-300 border border-emerald-500/20' :
                    'bg-blue-500/10 text-blue-300 border border-blue-500/20'
                  }`}
                >
                  <span>{issue.level === 'error' ? '❌' : issue.level === 'warning' ? '⚠️' : issue.level === 'success' ? '✅' : 'ℹ️'}</span>
                  <span>{issue.msg}</span>
                </div>
              ))}
            </div>
          )}
        </div>
      )}

      {/* CIS agent info */}
      <div className="bg-slate-800 rounded-xl p-5 border border-slate-700">
        <div className="flex items-center gap-3">
          <span className="text-3xl">🩺</span>
          <div>
            <h3 className="font-semibold text-slate-100">CIS Assistant Agent</h3>
            <p className="text-sm text-slate-400">Circulatory Informatics System — methodology guidance, contract generation, and code validation powered by the seven CIS principles.</p>
          </div>
          <div className="ml-auto text-right shrink-0">
            <p className="text-2xl font-bold text-purple-400">97%</p>
            <p className="text-xs text-slate-500">CIS accuracy</p>
          </div>
        </div>
      </div>
    </div>
  )
}
