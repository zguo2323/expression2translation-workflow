"""Filesystem-scoped project registry for E2T laboratory deployments."""
from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sqlite3


PROJECT_ID = re.compile(r"[a-z0-9][a-z0-9-]{0,62}\Z")
PROJECT_DIRECTORIES = ("inputs", "raw", "work", "results", "logs", "provenance")
STATUSES = {"created", "running", "succeeded", "failed", "archived"}
STAGES = {"init", "metadata", "qc", "rnaseq", "riboseq", "integration"}


class ProjectError(ValueError):
    """Raised for invalid project-registry operations."""


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def lab_root(path=None):
    """Return a normalized lab root without touching the filesystem."""
    return Path(path or "e2t-lab-data").expanduser().resolve()


def validate_project_id(project_id):
    if not isinstance(project_id, str) or not PROJECT_ID.fullmatch(project_id):
        raise ProjectError("Project ID must use 1-63 lowercase letters, digits, or hyphens and start with a letter or digit")
    return project_id


def registry_path(root):
    return lab_root(root) / "registry.sqlite"


SCHEMA = """
PRAGMA foreign_keys = ON;
CREATE TABLE IF NOT EXISTS projects (
  project_id TEXT PRIMARY KEY,
  title TEXT NOT NULL,
  owner TEXT,
  status TEXT NOT NULL CHECK(status IN ('created','running','succeeded','failed','archived')),
  project_path TEXT NOT NULL UNIQUE,
  created_at TEXT NOT NULL,
  updated_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS status_history (
  id INTEGER PRIMARY KEY,
  project_id TEXT NOT NULL,
  stage TEXT NOT NULL,
  status TEXT NOT NULL CHECK(status IN ('created','running','succeeded','failed','archived')),
  message TEXT NOT NULL,
  recorded_at TEXT NOT NULL,
  FOREIGN KEY(project_id) REFERENCES projects(project_id) ON DELETE CASCADE
);
CREATE TABLE IF NOT EXISTS artifacts (
  project_id TEXT NOT NULL,
  relative_path TEXT NOT NULL,
  kind TEXT NOT NULL,
  bytes INTEGER NOT NULL CHECK(bytes >= 0),
  sha256 TEXT NOT NULL CHECK(length(sha256) = 64),
  recorded_at TEXT NOT NULL,
  PRIMARY KEY(project_id, relative_path),
  FOREIGN KEY(project_id) REFERENCES projects(project_id) ON DELETE CASCADE
);
CREATE INDEX IF NOT EXISTS idx_status_history_project ON status_history(project_id, id);
"""


@contextmanager
def connection(root):
    root = lab_root(root)
    root.mkdir(parents=True, exist_ok=True)
    database = registry_path(root)
    conn = sqlite3.connect(database)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    conn.executescript(SCHEMA)
    try:
        yield conn
    finally:
        conn.close()


def project_path(root, project_id):
    project_id = validate_project_id(project_id)
    return lab_root(root) / "projects" / project_id


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def init_project(root, project_id, title=None, owner=None):
    """Create an isolated project directory and its registry record."""
    project_id = validate_project_id(project_id)
    root, target = lab_root(root), project_path(root, project_id)
    with connection(root) as conn:
        existing = conn.execute("SELECT 1 FROM projects WHERE project_id=?", (project_id,)).fetchone()
        if existing or target.exists():
            raise ProjectError(f"Project already exists: {project_id}")
        target.mkdir(parents=True)
        for name in PROJECT_DIRECTORIES:
            (target / name).mkdir()
        created = timestamp()
        metadata = {
            "schema_version": 1,
            "project_id": project_id,
            "title": title or project_id,
            "owner": owner,
            "created_at": created,
            "storage": {"raw": "raw", "work": "work", "results": "results", "logs": "logs", "provenance": "provenance"},
        }
        (target / "project.json").write_text(json.dumps(metadata, indent=2) + "\n")
        with conn:
            conn.execute("INSERT INTO projects VALUES (?,?,?,?,?,?,?)",
                         (project_id, metadata["title"], owner, "created", str(target), created, created))
            conn.execute("INSERT INTO status_history(project_id,stage,status,message,recorded_at) VALUES (?,?,?,?,?)",
                         (project_id, "init", "created", "Project initialized", created))
    return project_details(root, project_id)


def _project(conn, project_id):
    project_id = validate_project_id(project_id)
    row = conn.execute("SELECT * FROM projects WHERE project_id=?", (project_id,)).fetchone()
    if row is None:
        raise ProjectError(f"Unknown project: {project_id}")
    return row


def list_projects(root):
    with connection(root) as conn:
        return [dict(row) for row in conn.execute("SELECT * FROM projects ORDER BY created_at, project_id")]


def project_details(root, project_id):
    with connection(root) as conn:
        project = dict(_project(conn, project_id))
        project["history"] = [dict(row) for row in conn.execute(
            "SELECT stage,status,message,recorded_at FROM status_history WHERE project_id=? ORDER BY id", (project_id,))]
        project["artifacts"] = [dict(row) for row in conn.execute(
            "SELECT relative_path,kind,bytes,sha256,recorded_at FROM artifacts WHERE project_id=? ORDER BY relative_path", (project_id,))]
        return project


def set_status(root, project_id, stage, status, message=""):
    if stage not in STAGES:
        raise ProjectError(f"Unsupported stage: {stage}")
    if status not in STATUSES:
        raise ProjectError(f"Unsupported status: {status}")
    with connection(root) as conn, conn:
        _project(conn, project_id)
        recorded = timestamp()
        conn.execute("UPDATE projects SET status=?,updated_at=? WHERE project_id=?", (status, recorded, project_id))
        conn.execute("INSERT INTO status_history(project_id,stage,status,message,recorded_at) VALUES (?,?,?,?,?)",
                     (project_id, stage, status, message or "", recorded))
    return project_details(root, project_id)


def add_artifact(root, project_id, path, kind):
    """Register an existing file located inside the project's own directory."""
    if not isinstance(kind, str) or not kind.strip():
        raise ProjectError("Artifact kind must be non-empty")
    with connection(root) as conn, conn:
        project = _project(conn, project_id)
        base = Path(project["project_path"]).resolve()
        artifact = Path(path).expanduser().resolve()
        if not artifact.is_file():
            raise ProjectError(f"Artifact is not a file: {artifact}")
        try:
            relative = artifact.relative_to(base)
        except ValueError as error:
            raise ProjectError("Artifact must be inside its project directory") from error
        recorded = timestamp()
        conn.execute("INSERT INTO artifacts VALUES (?,?,?,?,?,?) "
                     "ON CONFLICT(project_id,relative_path) DO UPDATE SET kind=excluded.kind,bytes=excluded.bytes,sha256=excluded.sha256,recorded_at=excluded.recorded_at",
                     (project_id, str(relative), kind.strip(), artifact.stat().st_size, sha256(artifact), recorded))
    return project_details(root, project_id)


def check_registry(root):
    with connection(root) as conn:
        foreign = conn.execute("PRAGMA foreign_key_check").fetchall()
        quick = conn.execute("PRAGMA quick_check").fetchone()[0]
    if foreign or quick != "ok":
        raise ProjectError(f"Registry integrity failed: {foreign or quick}")
    return "ok"
