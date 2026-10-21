#!/usr/bin/env bash
# Runs the check of the note. Needs Python 3.9 or later and nothing else.
set -euo pipefail
cd "$(dirname "$0")"
python3 -I mirror_clause.py
