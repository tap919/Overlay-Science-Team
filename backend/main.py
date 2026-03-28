#!/usr/bin/env python3
"""
Overlay Science Team - Backend API
Full-stack agentic science team with CIS Assistant integration
"""
import asyncio
import importlib.util
import json
import os
import sys
import uuid
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Any
import logging

from fastapi import FastAPI, WebSocket, UploadFile, File, HTTPException, WebSocketDisconnect
from fastapi.middleware.cors import CORSMiddleware
import uvicorn

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

app = FastAPI(
    title="Overlay Science Team API",
    description="Full-stack agentic cross-disciplinary science team with CIS Assistant, cross-disciplinary agents, and circulatory informatics",
    version="1.0.0"
)

# CORS configuration: use an explicit allowlist, configurable via environment variables.
cors_origins_env = os.getenv("CORS_ALLOW_ORIGINS")
if cors_origins_env:
    allowed_origins = [origin.strip() for origin in cors_origins_env.split(",") if origin.strip()]
else:
    # Safe defaults for local development; adjust via CORS_ALLOW_ORIGINS in production.
    allowed_origins = ["http://localhost", "http://localhost:3000"]

allow_credentials = os.getenv("CORS_ALLOW_CREDENTIALS", "false").lower() == "true"

app.add_middleware(
    CORSMiddleware,
    allow_origins=allowed_origins,
    allow_credentials=allow_credentials,
    allow_methods=["*"],
    allow_headers=["*"],
)

WORKSPACE = Path("./science_workspace")
WORKSPACE.mkdir(exist_ok=True)
REPO_ROOT = Path(__file__).resolve().parent.parent

# In-memory storage
studies_db: Dict[str, Any] = {}
executions_db: Dict[str, Any] = {}
active_connections: Dict[str, List[WebSocket]] = {}

# Lock for safe concurrent updates to shared agent state
_agents_lock = asyncio.Lock()

# --- Science Team Agents ---
AGENTS = [
    {"id": "source-ingester", "name": "Source Ingester", "role": "Parses and extracts data from uploaded research sources (PDFs, CSVs, images)", "skills": {"parsing": 0.9, "extraction": 0.85, "ocr": 0.7}, "status": "idle", "memory": {"learned_workflows": 12, "successful_tasks": 48, "avg_execution_time_ms": 1200}},
    {"id": "quantum-runner", "name": "Quantum Chaos Runner", "role": "Runs quantum chaos simulations via GLCCE for protein folding and molecular analysis", "skills": {"quantum_simulation": 0.88, "chaos_analysis": 0.82, "optimization": 0.79}, "status": "idle", "memory": {"learned_workflows": 8, "successful_tasks": 31, "avg_execution_time_ms": 8500}},
    {"id": "biotech-simulator", "name": "Biotech Simulator", "role": "Executes molecular dynamics in Digital Lab, validates against biotech databases", "skills": {"molecular_dynamics": 0.91, "pathway_analysis": 0.87, "drug_binding": 0.83}, "status": "idle", "memory": {"learned_workflows": 15, "successful_tasks": 67, "avg_execution_time_ms": 12000}},
    {"id": "analytics-compiler", "name": "Analytics Compiler", "role": "Synthesizes insights from multi-source data using DataLite analytics engine", "skills": {"data_synthesis": 0.94, "statistical_analysis": 0.89, "visualization": 0.85}, "status": "idle", "memory": {"learned_workflows": 20, "successful_tasks": 89, "avg_execution_time_ms": 3200}},
    {"id": "writer-agent", "name": "Scientific Writer", "role": "Generates research papers and book chapters via Book Lab with structured scientific reasoning", "skills": {"paper_writing": 0.92, "hypothesis_generation": 0.88, "citation_management": 0.84}, "status": "idle", "memory": {"learned_workflows": 18, "successful_tasks": 72, "avg_execution_time_ms": 6700}},
    {"id": "cis-agent", "name": "CIS Assistant", "role": "Circulatory Informatics System - contract generation, validation, and methodology guidance", "skills": {"contract_generation": 0.95, "code_validation": 0.93, "cis_methodology": 0.97}, "status": "idle", "memory": {"learned_workflows": 25, "successful_tasks": 112, "avg_execution_time_ms": 800}},
    {"id": "metaphor-agent", "name": "Metaphor Engine", "role": "Enhances scientific content with comic metaphors and visual storytelling for accessibility", "skills": {"metaphor_generation": 0.87, "visual_storytelling": 0.82, "accessibility": 0.91}, "status": "idle", "memory": {"learned_workflows": 10, "successful_tasks": 44, "avg_execution_time_ms": 2100}},
]

# Cross-disciplinary disciplines
DISCIPLINES = [
    {"id": "quantum-physics", "name": "Quantum Physics", "color": "#6366f1", "agents": ["quantum-runner"], "active": True},
    {"id": "molecular-biology", "name": "Molecular Biology", "color": "#10b981", "agents": ["biotech-simulator"], "active": True},
    {"id": "data-science", "name": "Data Science", "color": "#f59e0b", "agents": ["analytics-compiler"], "active": True},
    {"id": "computational-chemistry", "name": "Computational Chemistry", "color": "#ef4444", "agents": ["biotech-simulator", "quantum-runner"], "active": True},
    {"id": "informatics", "name": "Circulatory Informatics", "color": "#8b5cf6", "agents": ["cis-agent"], "active": True},
    {"id": "science-writing", "name": "Science Communication", "color": "#06b6d4", "agents": ["writer-agent", "metaphor-agent"], "active": True},
]

PIPELINE_STAGES = [
    {"id": "ingest", "name": "Source Ingestion", "agent_id": "source-ingester", "order": 1},
    {"id": "quantum", "name": "Quantum Chaos Analysis", "agent_id": "quantum-runner", "order": 2},
    {"id": "biotech", "name": "Biotech Simulation", "agent_id": "biotech-simulator", "order": 3},
    {"id": "analytics", "name": "Insight Synthesis", "agent_id": "analytics-compiler", "order": 4},
    {"id": "writing", "name": "Paper & Book Writing", "agent_id": "writer-agent", "order": 5},
    {"id": "metaphor", "name": "Metaphor Enhancement", "agent_id": "metaphor-agent", "order": 6},
]


@lru_cache(maxsize=1)
def _load_digital_lab_tools() -> List[dict]:
    tools_path = REPO_ROOT / "Digital Lab tools"
    try:
        with tools_path.open("r", encoding="utf-8") as handle:
            return json.load(handle)
    except (FileNotFoundError, json.JSONDecodeError):
        logging.exception("Failed to load digital lab tools from %s", tools_path)
        return []


@lru_cache(maxsize=1)
def _load_enhanced_api_catalog() -> dict:
    module_path = REPO_ROOT / "enhanced-scientific-apis.py"
    try:
        spec = importlib.util.spec_from_file_location("overlay_enhanced_scientific_apis", module_path)
        if spec is None or spec.loader is None:
            raise ImportError(f"Unable to load scientific API module from {module_path}")

        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)

        catalog = module.get_enhanced_apis()
        apis = []
        category_counts: Dict[str, int] = {}

        for api_id, api in sorted(catalog.apis.items(), key=lambda item: item[1].name.lower()):
            category = api.category.value
            category_counts[category] = category_counts.get(category, 0) + 1
            apis.append({
                "id": api_id,
                "name": api.name,
                "category": category,
                "description": api.description,
                "install": api.install_command,
                "example": api.example_code,
                "docs": api.documentation_url,
                "github": api.github_url,
                "requires_api_key": api.requires_api_key,
            })

        categories = [
            {"id": category, "name": category.replace("_", " ").title(), "count": count}
            for category, count in sorted(category_counts.items())
        ]
        return {"apis": apis, "categories": categories}
    except Exception:
        logging.exception("Failed to load enhanced scientific APIs from %s", module_path)
        return {"apis": [], "categories": []}


@app.get("/")
async def root():
    return {"service": "Overlay Science Team API", "version": "1.0.0", "status": "online"}


@app.get("/api/v1/agents")
async def get_agents():
    return {"agents": AGENTS, "total": len(AGENTS)}


@app.get("/api/v1/disciplines")
async def get_disciplines():
    return {"disciplines": DISCIPLINES, "total": len(DISCIPLINES)}


@app.get("/api/v1/pipeline/stages")
async def get_pipeline_stages():
    return {"stages": PIPELINE_STAGES}


@app.post("/api/v1/upload")
async def upload_files(study_id: str, files: List[UploadFile] = File(...)):
    study_dir = WORKSPACE / study_id
    study_dir.mkdir(parents=True, exist_ok=True)

    uploaded = []
    for file in files:
        # Sanitize filename to prevent path traversal
        original_filename = file.filename or ""
        safe_filename = Path(original_filename).name
        if not safe_filename or safe_filename in (".", ".."):
            raise HTTPException(status_code=400, detail="Invalid filename")

        content = await file.read()
        file_path = study_dir / safe_filename
        with open(file_path, "wb") as f:
            f.write(content)
        uploaded.append({"filename": safe_filename, "size": len(content), "path": str(file_path)})

    studies_db[study_id] = {
        "id": study_id,
        "created_at": datetime.now().isoformat(),
        "files": uploaded,
        "status": "ready"
    }

    return {"study_id": study_id, "uploaded_files": uploaded, "status": "ready"}


@app.post("/api/v1/execute")
async def execute_pipeline(body: dict):
    study_id = body.get("study_id", f"study_{uuid.uuid4().hex[:8]}")
    study_type = body.get("study_type", "quantum_biotech")
    config = body.get("config", {})

    execution_id = str(uuid.uuid4())
    executions_db[execution_id] = {
        "id": execution_id,
        "study_id": study_id,
        "study_type": study_type,
        "config": config,
        "status": "queued",
        "progress": 0,
        "current_stage": "queued",
        "message": "Pipeline queued",
        "created_at": datetime.now().isoformat(),
        "stages": {}
    }

    # Run pipeline in background
    asyncio.create_task(run_pipeline(execution_id, study_id, study_type))

    return {"execution_id": execution_id, "study_id": study_id, "status": "queued"}


@app.get("/api/v1/executions/{execution_id}")
async def get_execution(execution_id: str):
    if execution_id not in executions_db:
        raise HTTPException(status_code=404, detail="Execution not found")
    return executions_db[execution_id]


@app.get("/api/v1/studies")
async def get_studies():
    return {"studies": list(studies_db.values()), "total": len(studies_db)}


@app.websocket("/ws/executions/{execution_id}")
async def websocket_endpoint(websocket: WebSocket, execution_id: str):
    await websocket.accept()

    # Validate that the execution exists before registering the connection
    if execution_id not in executions_db:
        await websocket.send_json({"error": "Execution not found", "execution_id": execution_id})
        await websocket.close()
        return

    if execution_id not in active_connections:
        active_connections[execution_id] = []
    active_connections[execution_id].append(websocket)

    try:
        # Stream execution updates until a terminal status is reached or the execution disappears
        while True:
            exec_data = executions_db.get(execution_id)
            if not exec_data:
                # Execution no longer tracked; stop streaming
                break

            await websocket.send_json(exec_data)

            status = exec_data.get("status")
            if status in {"success", "failed", "error", "cancelled"}:
                # Execution has reached a terminal state; close the connection
                break

            await asyncio.sleep(1)
    except WebSocketDisconnect:
        # Client disconnected; cleanup handled in finally
        pass
    finally:
        # Ensure the websocket is removed from active connections
        connections = active_connections.get(execution_id)
        if connections and websocket in connections:
            connections.remove(websocket)
            if not connections:
                active_connections.pop(execution_id, None)
        try:
            await websocket.close()
        except Exception:
            # Ignore errors during cleanup close
            pass


async def run_pipeline(execution_id: str, study_id: str, study_type: str):
    """Simulate pipeline execution through all stages"""
    exec_data = executions_db[execution_id]
    exec_data["status"] = "running"
    exec_data["started_at"] = datetime.now().isoformat()

    for i, stage in enumerate(PIPELINE_STAGES):
        stage_id = stage["id"]
        exec_data["current_stage"] = stage["name"]
        exec_data["message"] = f"Running {stage['name']}..."
        exec_data["stages"][stage_id] = {"status": "running", "progress": 0, "started_at": datetime.now().isoformat()}

        # Update agent status — protected by lock to prevent concurrent-execution races
        async with _agents_lock:
            for agent in AGENTS:
                if agent["id"] == stage["agent_id"]:
                    agent["status"] = "busy"

        # Simulate stage execution in steps
        for pct in range(0, 101, 20):
            await asyncio.sleep(0.5)
            exec_data["stages"][stage_id]["progress"] = pct
            exec_data["progress"] = int(((i + pct / 100) / len(PIPELINE_STAGES)) * 100)

        exec_data["stages"][stage_id]["status"] = "complete"
        exec_data["stages"][stage_id]["progress"] = 100
        exec_data["stages"][stage_id]["completed_at"] = datetime.now().isoformat()
        exec_data["stages"][stage_id]["result"] = _generate_stage_result(stage_id, study_id)

        # Reset agent status — protected by lock
        async with _agents_lock:
            for agent in AGENTS:
                if agent["id"] == stage["agent_id"]:
                    agent["status"] = "idle"
                    agent["memory"]["successful_tasks"] += 1

    exec_data["status"] = "success"
    exec_data["progress"] = 100
    exec_data["current_stage"] = "complete"
    exec_data["message"] = "Pipeline complete"
    exec_data["completed_at"] = datetime.now().isoformat()
    exec_data["outputs"] = _generate_outputs(study_id)


def _generate_stage_result(stage_id: str, study_id: str) -> dict:
    results = {
        "ingest": {"sources_parsed": 3, "entities_extracted": 142, "quality_score": 0.94},
        "quantum": {"convergence_rate": 0.87, "chaos_metric": 0.63, "folding_energy": -234.5, "stable_states": 12},
        "biotech": {"binding_affinity": 8.2, "pathway_activation": ["MAPK", "PI3K", "mTOR"], "toxicity_score": 0.12},
        "analytics": {"insights_generated": 18, "correlations_found": 7, "confidence_avg": 0.89},
        "writing": {"papers_drafted": 2, "book_chapters": 1, "total_pages": 44, "citations": 38},
        "metaphor": {"metaphors_added": 12, "accessibility_score": 0.91, "visual_aids": 8},
    }
    return results.get(stage_id, {})


def _generate_outputs(study_id: str) -> list:
    ts = datetime.now().isoformat()
    return [
        {"id": f"{study_id}_paper1", "type": "paper", "title": "Quantum Chaos Analysis of Protein Folding Dynamics", "timestamp": ts, "pages": 12, "citations": 23, "status": "ready"},
        {"id": f"{study_id}_paper2", "type": "paper", "title": "Biotech Simulation Results: Novel Pathways in Molecular Binding", "timestamp": ts, "pages": 8, "citations": 15, "status": "ready"},
        {"id": f"{study_id}_book", "type": "book_chapter", "title": "Chapter 3: Chaos Theory in Biological Systems", "timestamp": ts, "pages": 24, "metaphors": 8, "status": "ready"},
        {"id": f"{study_id}_analytics", "type": "analytics", "title": "Convergence Analysis Dashboard", "timestamp": ts, "charts": 12, "datasets": 4, "status": "ready"},
    ]


@app.get("/api/v1/cis/principles")
async def get_cis_principles():
    """CIS Seven Principles for circulatory informatics"""
    principles = [
        {"id": "distributed_autonomy", "name": "Distributed Autonomy", "description": "No single point of control. Each organ operates autonomously within a decentralized governance framework.", "icon": "🔄"},
        {"id": "continuous_sensing", "name": "Continuous Sensing", "description": "The organism must continuously monitor its own state through events. Observability is part of the nervous system.", "icon": "📡"},
        {"id": "feedback_driven_adaptation", "name": "Feedback-Driven Adaptation", "description": "Adaptation requires feedback loops that measure deviation and trigger corrective action.", "icon": "🔁"},
        {"id": "emergent_intelligence", "name": "Emergent Intelligence", "description": "Intelligence emerges from interaction of simple agents following local rules.", "icon": "🧠"},
        {"id": "memory_and_learning", "name": "Memory and Learning", "description": "The organism learns from past events and adjusts future behavior.", "icon": "💾"},
        {"id": "graceful_degradation", "name": "Graceful Degradation", "description": "No single component is essential. The system continues functioning when components fail.", "icon": "🛡️"},
        {"id": "efficient_resource_flow", "name": "Efficient Resource Flow", "description": "Resources flow to where they're needed. Demand drives allocation.", "icon": "⚡"},
    ]
    return {"principles": principles}


@app.get("/api/v1/cis/capabilities")
async def get_cis_capabilities():
    digital_lab_tools = _load_digital_lab_tools()
    api_catalog = _load_enhanced_api_catalog()
    return {
        "digital_lab_tools": digital_lab_tools,
        "enhanced_apis": api_catalog["apis"],
        "api_categories": api_catalog["categories"],
        "summary": {
            "digital_lab_tool_count": len(digital_lab_tools),
            "enhanced_api_count": len(api_catalog["apis"]),
            "api_category_count": len(api_catalog["categories"]),
        }
    }


@app.get("/api/v1/health")
async def health_check():
    return {"status": "healthy", "timestamp": datetime.now().isoformat(), "agents_active": sum(1 for a in AGENTS if a["status"] == "busy")}


if __name__ == "__main__":
    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True)
