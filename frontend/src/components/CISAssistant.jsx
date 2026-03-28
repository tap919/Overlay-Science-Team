import React, { useState, useEffect } from 'react'

export default function CISAssistant({ apiBase }) {
  const [principles, setPrinciples] = useState([])
  const [capabilities, setCapabilities] = useState({
    digital_lab_tools: [],
    enhanced_apis: [],
    api_categories: [],
    summary: null,
  })
  const [selectedPrinciple, setSelectedPrinciple] = useState(null)
  const [contractInput, setContractInput] = useState('')
  const [contractOutput, setContractOutput] = useState('')
  const [validationInput, setValidationInput] = useState('')
  const [validationResult, setValidationResult] = useState(null)
  const [activeSection, setActiveSection] = useState('principles')

  useEffect(() => {
    Promise.all([
      fetch(`${apiBase}/cis/principles`).then(r => r.json()),
      fetch(`${apiBase}/cis/capabilities`).then(r => r.json()),
    ])
      .then(([principlesData, capabilitiesData]) => {
        setPrinciples(principlesData.principles || [])
        setCapabilities({
          digital_lab_tools: capabilitiesData.digital_lab_tools || [],
          enhanced_apis: capabilitiesData.enhanced_apis || [],
          api_categories: capabilitiesData.api_categories || [],
          summary: capabilitiesData.summary || null,
        })
      })
      .catch(console.error)
  }, [apiBase])

  const generateContract = () => {
    if (!contractInput.trim()) return
    // Sanitize to a valid TypeScript identifier:
    // 1. Replace whitespace and hyphens with underscores
    // 2. Strip any remaining non-alphanumeric/underscore characters
    // 3. Prepend underscore if the result starts with a digit
    const rawName = contractInput.trim().replace(/[\s-]+/g, '_').replace(/[^a-zA-Z0-9_]/g, '')
    const safeName = /^\d/.test(rawName) ? `_${rawName}` : (rawName || 'Unnamed')
    // Simulate contract generation
    setContractOutput(`/**
 * CIS Contract: ${contractInput.trim()}
 * Generated: ${new Date().toISOString()}
 * Principle: Distributed Autonomy + Feedback-Driven Adaptation
 */
interface ${safeName}Contract {
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
    { id: 'tooling', label: '🧪 Lab Tools & APIs' },
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
      {activeSection === 'tooling' && (
        <div className="space-y-6">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            {[
              {
                label: 'Digital Lab Tools',
                value: capabilities.summary?.digital_lab_tool_count ?? capabilities.digital_lab_tools.length,
                tone: 'text-emerald-300',
              },
              {
                label: 'Enhanced APIs',
                value: capabilities.summary?.enhanced_api_count ?? capabilities.enhanced_apis.length,
                tone: 'text-indigo-300',
              },
              {
                label: 'API Categories',
                value: capabilities.summary?.api_category_count ?? capabilities.api_categories.length,
                tone: 'text-amber-300',
              },
            ].map(card => (
              <div key={card.label} className="bg-slate-800 rounded-xl p-5 border border-slate-700">
                <p className="text-sm text-slate-400">{card.label}</p>
                <p className={`text-3xl font-bold mt-2 ${card.tone}`}>{card.value}</p>
              </div>
            ))}
          </div>

          <div className="bg-slate-800 rounded-xl p-6 border border-slate-700">
            <h2 className="text-lg font-semibold text-indigo-300">🧬 Digital Lab Tooling</h2>
            <p className="text-sm text-slate-400 mt-1">
              Curated biotech and bioinformatics tools ready to support CIS-guided digital lab workflows.
            </p>
            <div className="mt-4 grid grid-cols-1 md:grid-cols-2 gap-4">
              {capabilities.digital_lab_tools.map(tool => (
                <div key={tool.name} className="rounded-xl border border-slate-700 bg-slate-900/60 p-4 space-y-3">
                  <div>
                    <div className="flex items-center justify-between gap-3">
                      <h3 className="font-semibold text-slate-100">{tool.name}</h3>
                      <span className="text-xs px-2 py-1 rounded-full bg-emerald-500/10 text-emerald-300 border border-emerald-500/20">
                        {tool.category}
                      </span>
                    </div>
                    <p className="text-sm text-slate-400 mt-2">{tool.description}</p>
                  </div>
                  <div className="text-xs text-slate-500">
                    <span className="text-slate-400">Install:</span> {tool.install}
                  </div>
                  <pre className="bg-slate-950 rounded-lg p-3 text-xs text-emerald-300 overflow-x-auto whitespace-pre-wrap border border-slate-800">
                    {tool.example}
                  </pre>
                  <div className="flex gap-4 text-sm">
                    <a className="text-indigo-300 hover:text-indigo-200" href={tool.docs} target="_blank" rel="noreferrer">
                      Documentation ↗
                    </a>
                    <a className="text-slate-300 hover:text-slate-200" href={tool.github} target="_blank" rel="noreferrer">
                      GitHub ↗
                    </a>
                  </div>
                </div>
              ))}
            </div>
          </div>

          <div className="bg-slate-800 rounded-xl p-6 border border-slate-700">
            <div className="flex flex-wrap items-center gap-2">
              <h2 className="text-lg font-semibold text-indigo-300 mr-2">🔌 Enhanced Scientific APIs</h2>
              {capabilities.api_categories.map(category => (
                <span
                  key={category.id}
                  className="text-xs px-2 py-1 rounded-full bg-indigo-500/10 text-indigo-300 border border-indigo-500/20"
                >
                  {category.name} · {category.count}
                </span>
              ))}
            </div>
            <p className="text-sm text-slate-400 mt-3">
              Expanded API coverage for compounds, literature, genes, proteins, chemistry, and genomics workflows.
            </p>
            <div className="mt-4 space-y-3">
              {capabilities.enhanced_apis.map(api => (
                <details
                  key={api.id}
                  className="rounded-xl border border-slate-700 bg-slate-900/60 p-4 group"
                >
                  <summary className="cursor-pointer list-none flex flex-col md:flex-row md:items-center md:justify-between gap-2">
                    <div>
                      <h3 className="font-semibold text-slate-100">{api.name}</h3>
                      <p className="text-sm text-slate-400 mt-1">{api.description}</p>
                    </div>
                    <div className="flex items-center gap-2 text-xs">
                      <span className="px-2 py-1 rounded-full bg-blue-500/10 text-blue-300 border border-blue-500/20">
                        {api.category.replace('_', ' ')}
                      </span>
                      {api.requires_api_key && (
                        <span className="px-2 py-1 rounded-full bg-amber-500/10 text-amber-300 border border-amber-500/20">
                          API key required
                        </span>
                      )}
                    </div>
                  </summary>
                  <div className="mt-4 space-y-3">
                    <div className="text-xs text-slate-500">
                      <span className="text-slate-400">Install:</span> {api.install}
                    </div>
                    <pre className="bg-slate-950 rounded-lg p-3 text-xs text-emerald-300 overflow-x-auto whitespace-pre-wrap border border-slate-800">
                      {api.example}
                    </pre>
                    <div className="flex gap-4 text-sm">
                      <a className="text-indigo-300 hover:text-indigo-200" href={api.docs} target="_blank" rel="noreferrer">
                        Documentation ↗
                      </a>
                      <a className="text-slate-300 hover:text-slate-200" href={api.github} target="_blank" rel="noreferrer">
                        GitHub ↗
                      </a>
                    </div>
                  </div>
                </details>
              ))}
            </div>
          </div>
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
