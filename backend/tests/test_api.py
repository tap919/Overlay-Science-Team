"""Tests for the Overlay Science Team API"""
from pathlib import Path
from fastapi.testclient import TestClient
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from main import app

client = TestClient(app)


def test_root():
    resp = client.get("/")
    assert resp.status_code == 200
    assert resp.json()["service"] == "Overlay Science Team API"


def test_get_agents():
    resp = client.get("/api/v1/agents")
    assert resp.status_code == 200
    data = resp.json()
    assert "agents" in data
    assert len(data["agents"]) >= 6


def test_get_disciplines():
    resp = client.get("/api/v1/disciplines")
    assert resp.status_code == 200
    data = resp.json()
    assert "disciplines" in data


def test_get_pipeline_stages():
    resp = client.get("/api/v1/pipeline/stages")
    assert resp.status_code == 200
    data = resp.json()
    assert "stages" in data
    assert len(data["stages"]) == 6


def test_get_cis_principles():
    resp = client.get("/api/v1/cis/principles")
    assert resp.status_code == 200
    data = resp.json()
    assert "principles" in data
    assert len(data["principles"]) == 7


def test_health():
    resp = client.get("/api/v1/health")
    assert resp.status_code == 200
    assert resp.json()["status"] == "healthy"


def test_execute_pipeline():
    resp = client.post("/api/v1/execute", json={"study_id": "test-study-001", "study_type": "quantum_biotech"})
    assert resp.status_code == 200
    data = resp.json()
    assert "execution_id" in data
    assert data["status"] == "queued"
