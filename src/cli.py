"""Terminal entry point for E2T laboratory project management."""
from __future__ import annotations

import argparse
import json
import os
import sys

from src.lab.registry import ProjectError, add_artifact, check_registry, init_project, list_projects, project_details, set_status


def root_argument(parser):
    parser.add_argument("--root", default=os.environ.get("E2T_HOME", "e2t-lab-data"),
                        help="Laboratory storage root (default: $E2T_HOME or ./e2t-lab-data)")


def parser():
    command = argparse.ArgumentParser(prog="e2t", description=__doc__)
    subcommands = command.add_subparsers(dest="command", required=True)

    init = subcommands.add_parser("init", help="Create an isolated project")
    init.add_argument("project_id")
    init.add_argument("--title")
    init.add_argument("--owner")
    root_argument(init)

    listing = subcommands.add_parser("list", help="List registered projects")
    listing.add_argument("--json", action="store_true")
    root_argument(listing)

    status = subcommands.add_parser("status", help="Show project history and registered artifacts")
    status.add_argument("project_id")
    status.add_argument("--json", action="store_true")
    root_argument(status)

    update = subcommands.add_parser("set-status", help="Record a project stage and status")
    update.add_argument("project_id")
    update.add_argument("--stage", required=True, choices=("init", "metadata", "qc", "rnaseq", "riboseq", "integration"))
    update.add_argument("--status", required=True, choices=("created", "running", "succeeded", "failed", "archived"))
    update.add_argument("--message", default="")
    root_argument(update)

    artifact = subcommands.add_parser("add-artifact", help="Register a result file inside one project")
    artifact.add_argument("project_id")
    artifact.add_argument("--path", required=True)
    artifact.add_argument("--kind", required=True)
    root_argument(artifact)

    check = subcommands.add_parser("check", help="Check registry SQLite integrity")
    root_argument(check)
    return command


def emit(value, as_json=False):
    if as_json:
        print(json.dumps(value, ensure_ascii=False, indent=2))
        return
    if isinstance(value, list):
        if not value:
            print("No registered projects.")
            return
        print("project_id\tstatus\towner\tproject_path")
        for row in value:
            print(f"{row['project_id']}\t{row['status']}\t{row['owner'] or '—'}\t{row['project_path']}")
        return
    print(f"{value['project_id']}: {value['status']} ({value['project_path']})")


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        if args.command == "init":
            result = init_project(args.root, args.project_id, args.title, args.owner)
            emit(result)
        elif args.command == "list":
            emit(list_projects(args.root), args.json)
        elif args.command == "status":
            emit(project_details(args.root, args.project_id), args.json)
        elif args.command == "set-status":
            emit(set_status(args.root, args.project_id, args.stage, args.status, args.message))
        elif args.command == "add-artifact":
            emit(add_artifact(args.root, args.project_id, args.path, args.kind))
        else:
            print(check_registry(args.root))
    except ProjectError as error:
        parser().error(str(error))
    return 0


if __name__ == "__main__":
    sys.exit(main())
