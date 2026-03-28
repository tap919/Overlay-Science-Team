import React, { useState } from 'react'

const SKILL_BAR_COLORS = ['bg-indigo-500', 'bg-emerald-500', 'bg-amber-500', 'bg-rose-500', 'bg-cyan-500']

function AgentCard({ agent, color }) {
  const skills = Object.entries(agent.skills || {})
  return (
    <div className="bg-slate-700 rounded-xl p-5 border border-slate-600 hover:border-slate-500 transition-colors">
      <div className="flex items-start justify-between mb-3">
        <div>
          <h3 className="font-semibold text-slate-100">{agent.name}</h3>
          <p className="text-xs text-slate-400 mt-0.5 leading-tight">{agent.role}</p>
        </div>
        <span
          className={`text-xs px-2 py-0.5 rounded-full font-medium shrink-0 ml-2 ${
            agent.status === 'busy' ? 'bg-blue-500/20 text-blue-300' : 'bg-emerald-500/20 text-emerald-300'
          }`}
        >
          {agent.status}
        </span>
      </div>

      {/* Skills */}
      <div className="space-y-1.5 mb-4">
        {skills.slice(0, 3).map(([skill, val], i) => (
          <div key={skill} className="flex items-center gap-2">
            <span className="text-xs text-slate-400 w-28 truncate capitalize">{skill.replace(/_/g, ' ')}</span>
            <div className="flex-1 bg-slate-600 rounded-full h-1.5">
              <div
                className={`h-1.5 rounded-full ${SKILL_BAR_COLORS[i % SKILL_BAR_COLORS.length]}`}
                style={{ width: `${Math.round(val * 100)}%` }}
              />
            </div>
            <span className="text-xs text-slate-400 w-8 text-right">{Math.round(val * 100)}%</span>
          </div>
        ))}
      </div>

      {/* Memory stats */}
      {agent.memory && (
        <div className="grid grid-cols-3 gap-2 text-center border-t border-slate-600 pt-3">
          <div>
            <p className="text-lg font-bold" style={{ color }}>{agent.memory.learned_workflows}</p>
            <p className="text-xs text-slate-500">workflows</p>
          </div>
          <div>
            <p className="text-lg font-bold" style={{ color }}>{agent.memory.successful_tasks}</p>
            <p className="text-xs text-slate-500">tasks done</p>
          </div>
          <div>
            <p className="text-lg font-bold" style={{ color }}>{(agent.memory.avg_execution_time_ms / 1000).toFixed(1)}s</p>
            <p className="text-xs text-slate-500">avg time</p>
          </div>
        </div>
      )}
    </div>
  )
}

export default function CrossDisciplinaryAgents({ agents, disciplines }) {
  const [selectedDiscipline, setSelectedDiscipline] = useState('all')

  const agentsForDiscipline = () => {
    if (selectedDiscipline === 'all') return agents
    const disc = disciplines.find(d => d.id === selectedDiscipline)
    if (!disc) return agents
    return agents.filter(a => disc.agents.includes(a.id))
  }

  const disciplineColor = (discId) => {
    const disc = disciplines.find(d => d.id === discId)
    return disc?.color || '#6366f1'
  }

  const agentColor = (agentId) => {
    for (const disc of disciplines) {
      if (disc.agents.includes(agentId)) return disc.color
    }
    return '#6366f1'
  }

  return (
    <div className="space-y-6">
      {/* Discipline filter */}
      <div className="bg-slate-800 rounded-xl p-6 border border-slate-700">
        <h2 className="text-lg font-semibold mb-4 text-indigo-300">🔬 Disciplines</h2>
        <div className="flex flex-wrap gap-2">
          <button
            onClick={() => setSelectedDiscipline('all')}
            className={`px-4 py-2 rounded-full text-sm font-medium transition-colors ${
              selectedDiscipline === 'all'
                ? 'bg-indigo-600 text-white'
                : 'bg-slate-700 text-slate-300 hover:bg-slate-600'
            }`}
          >
            All Agents ({agents.length})
          </button>
          {disciplines.map(disc => (
            <button
              key={disc.id}
              onClick={() => setSelectedDiscipline(disc.id)}
              className={`px-4 py-2 rounded-full text-sm font-medium transition-colors border ${
                selectedDiscipline === disc.id ? 'text-white' : 'bg-slate-700 text-slate-300 hover:bg-slate-600'
              }`}
              style={
                selectedDiscipline === disc.id
                  ? { backgroundColor: disc.color, borderColor: disc.color }
                  : { borderColor: disc.color + '44' }
              }
            >
              {disc.name}
            </button>
          ))}
        </div>
      </div>

      {/* Collaboration graph hint */}
      {selectedDiscipline !== 'all' && (() => {
        const disc = disciplines.find(d => d.id === selectedDiscipline)
        return disc ? (
          <div
            className="rounded-xl p-4 border text-sm"
            style={{ borderColor: disc.color + '55', backgroundColor: disc.color + '11' }}
          >
            <span style={{ color: disc.color }}>● {disc.name}</span>
            <span className="text-slate-400 ml-2">— {disc.agents.length} agent(s) assigned to this discipline</span>
          </div>
        ) : null
      })()}

      {/* Agent cards */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-5">
        {agentsForDiscipline().map(agent => (
          <AgentCard key={agent.id} agent={agent} color={agentColor(agent.id)} />
        ))}
      </div>

      {/* Cross-discipline collaboration matrix */}
      <div className="bg-slate-800 rounded-xl p-6 border border-slate-700">
        <h2 className="text-lg font-semibold mb-4 text-indigo-300">🔗 Cross-Discipline Collaboration</h2>
        <div className="overflow-x-auto">
          <table className="w-full text-sm">
            <thead>
              <tr>
                <th className="text-left text-slate-400 font-medium pb-2 pr-4">Discipline</th>
                <th className="text-left text-slate-400 font-medium pb-2 pr-4">Agents</th>
                <th className="text-left text-slate-400 font-medium pb-2">Shared With</th>
              </tr>
            </thead>
            <tbody>
              {disciplines.map(disc => {
                const shared = disciplines.filter(
                  d => d.id !== disc.id && d.agents.some(a => disc.agents.includes(a))
                )
                return (
                  <tr key={disc.id} className="border-t border-slate-700">
                    <td className="py-2 pr-4">
                      <span className="font-medium" style={{ color: disc.color }}>{disc.name}</span>
                    </td>
                    <td className="py-2 pr-4 text-slate-400">{disc.agents.join(', ')}</td>
                    <td className="py-2">
                      {shared.length > 0
                        ? shared.map(s => (
                            <span key={s.id} className="mr-2 text-xs px-2 py-0.5 rounded-full" style={{ backgroundColor: s.color + '22', color: s.color }}>
                              {s.name}
                            </span>
                          ))
                        : <span className="text-slate-600 text-xs">—</span>
                      }
                    </td>
                  </tr>
                )
              })}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  )
}
