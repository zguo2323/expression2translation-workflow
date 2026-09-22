import json
from pathlib import Path
import tempfile
import unittest

from src.lab.registry import ProjectError, add_artifact, check_registry, init_project, list_projects, project_details, set_status


class LabRegistryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name) / "lab"

    def test_init_status_artifact_and_integrity(self):
        project = init_project(self.root, "aging-yeast-2026", title="Aging yeast", owner="researcher")
        folder = Path(project["project_path"])
        self.assertTrue(all((folder / name).is_dir() for name in ("inputs", "raw", "work", "results", "logs", "provenance")))
        metadata = json.loads((folder / "project.json").read_text())
        self.assertEqual(metadata["project_id"], "aging-yeast-2026")
        set_status(self.root, "aging-yeast-2026", "qc", "running", "FastQC queued")
        report = folder / "results" / "report.html"
        report.write_text("report")
        details = add_artifact(self.root, "aging-yeast-2026", report, "html_report")
        self.assertEqual(details["status"], "running")
        self.assertEqual(details["artifacts"][0]["relative_path"], "results/report.html")
        self.assertEqual(len(details["history"]), 2)
        self.assertEqual(check_registry(self.root), "ok")

    def test_rejects_invalid_id_duplicates_and_external_artifact(self):
        with self.assertRaisesRegex(ProjectError, "Project ID"):
            init_project(self.root, "../escape")
        init_project(self.root, "valid-project")
        with self.assertRaisesRegex(ProjectError, "already exists"):
            init_project(self.root, "valid-project")
        external = Path(self.temporary.name) / "external.txt"
        external.write_text("no")
        with self.assertRaisesRegex(ProjectError, "inside its project"):
            add_artifact(self.root, "valid-project", external, "text")

    def test_list_and_unknown_project(self):
        init_project(self.root, "b-project")
        init_project(self.root, "a-project")
        self.assertEqual([row["project_id"] for row in list_projects(self.root)], ["b-project", "a-project"])
        with self.assertRaisesRegex(ProjectError, "Unknown project"):
            project_details(self.root, "missing-project")


if __name__ == "__main__":
    unittest.main()
