import React, { useState, useEffect } from 'react'

const CIS_PRINCIPLES_SHORT = [
  { name: 'Distributed Autonomy', icon: '🔄', color: '#6366f1' },
  { name: 'Continuous Sensing', icon: '📡', color: '#10b981' },
  { name: 'Feedback-Driven Adaptation', icon: '🔁', color: '#f59e0b' },
  { name: 'Emergent Intelligence', icon: '🧠', color: '#ef4444' },
  { name: 'Memory & Learning', icon: '💾', color: '#8b5cf6' },
  { name: 'Graceful Degradation', icon: '🛡️', color: '#06b6d4' },
  { name: 'Efficient Resource Flow', icon: '⚡', color: '#ec4899' },
]

function MetricCard({ label, value, unit, color, icon }) {
  return (
    <div className="bg-slate-700 rounded-xl p-5 border border-slate-600">
      <div className="flex items-center gap-2 mb-2">
        <span className="text-xl">{icon}</span>
        <span className="text-sm text-slate-400">{label}</span>
      </div>
      <p className="text-3xl font-bold" style={{ color }}>
        {value}<span className="text-base font-normal text-slate-400 ml-1">{unit}</span>
      </p>
    </div>
  )
}

export default function CirculatoryInformatics({ agents, health }) {
  const [tick, setTick] = useState(0)

  // Live-update tick for simulated vitals
  useEffect(() => {
    const id = setInterval(() => setTick(t => t + 1), 2000)
    return () => clearInterval(id)
  }, [])

  const idleAgents = agents.filter(a => a.status === 'idle').length
  const busyAgents = agents.filter(a => a.status === 'busy').length
  const totalTasks = agents.reduce((s, a) => s + (a.memory?.successful_tasks || 0), 0)
  const avgConfidence = agents.length
    ? (agents.reduce((s, a) => {
        const skills = Object.values(a.skills || {})
        return s + (skills.length ? skills.reduce((x, v) => x + v, 0) / skills.length : 0)
      }, 0) / agents.length * 100).toFixed(1)
    : '—'

  // Simulated vitals that shift slightly each tick
  const heartRate = 72 + Math.sin(tick * 0.7) * 3
  const bloodFlow = 85 + Math.sin(tick * 1.1) * 5
  const oxygenation = 97 + Math.sin(tick * 0.4) * 1

  return (
    <div className="space-y-6">
      {/* System vitals */}
      <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
        <MetricCard label="System Pulse" value={heartRate.toFixed(0)} unit="bpm" color="#ef4444" icon="❤️" />
        <MetricCard label="Data Flow" value={bloodFlow.toFixed(0)} unit="%" color="#10b981" icon="🩸" />
        <MetricCard label="Agent Oxygenation" value={oxygenation.toFixed(1)} unit="%" color="#6366f1" icon="💨" />
        <MetricCard label="Tasks Completed" value={totalTasks} unit="total" color="#f59e0b" icon="✅" />
      </div>

      {/* Agent circulation status */}
      <div className="bg-slate-800 rounded-xl p-6 border border-slate-700">
        <h2 className="text-lg font-semibold mb-4 text-indigo-300">🫀 Agent Circulatory Status</h2>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          {agents.map(agent => {
            const skills = Object.values(agent.skills || {})
            const avgSkill = skills.length ? skills.reduce((s, v) => s + v, 0) / skills.length : 0
            return (
              <div key={agent.id} className="flex items-center gap-4 bg-slate-700 rounded-lg px-4 py-3">
                {/* Pulse indicator */}
                <div className="relative w-10 h-10 shrink-0 flex items-center justify-center">
                  <div
                    className={`w-8 h-8 rounded-full ${agent.status === 'busy' ? 'animate-ping opacity-30' : ''}`}
                    style={{ backgroundColor: agent.status === 'busy' ? '#6366f1' : '#10b981' }}
                  />
                  <div
                    className="absolute w-6 h-6 rounded-full"
                    style={{ backgroundColor: agent.status === 'busy' ? '#6366f1' : '#10b981' }}
                  />
                </div>
                <div className="flex-1 min-w-0">
                  <div className="flex items-center justify-between">
                    <span className="text-sm font-medium text-slate-200 truncate">{agent.name}</span>
                    <span className={`text-xs ml-2 shrink-0 ${agent.status === 'busy' ? 'text-blue-400' : 'text-emerald-400'}`}>
                      {agent.status}
                    </span>
                  </div>
                  <div className="mt-1 w-full bg-slate-600 rounded-full h-1.5">
                    <div
                      className="h-1.5 rounded-full bg-indigo-500 transition-all duration-1000"
                      style={{ width: `${Math.round(avgSkill * 100)}%` }}
                    />
                  </div>
                  <div className="flex justify-between mt-0.5">
                    <span className="text-xs text-slate-500">{agent.memory?.successful_tasks || 0} tasks done</span>
                    <span className="text-xs text-slate-500">{Math.round(avgSkill * 100)}% capacity</span>
                  </div>
                </div>
              </div>
            )
          })}
        </div>

        {/* Summary */}
        <div className="mt-4 flex items-center gap-6 text-sm text-slate-400 border-t border-slate-700 pt-4">
          <span>🟢 {idleAgents} idle</span>
          <span>🔵 {busyAgents} active</span>
          <span>
            🧠 Avg confidence:{' '}
            {typeof avgConfidence === 'number' ? `${avgConfidence}%` : avgConfidence}
          </span>
        </div>
      </div>

      {/* CIS Principles health check */}
      <div className="bg-slate-800 rounded-xl p-6 border border-slate-700">
        <h2 className="text-lg font-semibold mb-4 text-indigo-300">🔬 CIS Principle Health</h2>
        <div className="space-y-3">
          {CIS_PRINCIPLES_SHORT.map((p, i) => {
            // Simulated health score per principle, drifts with tick
            const score = Math.min(100, Math.max(75, 90 + Math.sin((tick + i) * 0.5) * 8))
            return (
              <div key={p.name} className="flex items-center gap-3">
                <span className="text-lg w-6">{p.icon}</span>
                <span className="text-sm text-slate-300 w-52 shrink-0">{p.name}</span>
                <div className="flex-1 bg-slate-700 rounded-full h-2">
                  <div
                    className="h-2 rounded-full transition-all duration-1000"
                    style={{ width: `${score.toFixed(0)}%`, backgroundColor: p.color }}
                  />
                </div>
                <span className="text-xs text-slate-400 w-10 text-right">{score.toFixed(0)}%</span>
              </div>
            )
          })}
        </div>
      </div>

      {/* System health */}
      {health && (
        <div className="bg-slate-800 rounded-xl p-5 border border-slate-700 flex items-center gap-4">
          <span className="text-3xl">🏥</span>
          <div>
            <p className="font-semibold text-slate-100">System Health: <span className="text-emerald-400 capitalize">{health.status}</span></p>
            <p className="text-sm text-slate-400">Last checked: {new Date(health.timestamp).toLocaleTimeString()}</p>
          </div>
          <div className="ml-auto text-right">
            <p className="text-2xl font-bold text-emerald-400">{health.agents_active}</p>
            <p className="text-xs text-slate-500">agents active</p>
          </div>
        </div>
      )}
    </div>
  )
}
