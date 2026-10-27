#!/usr/bin/env bash
# Runs the checks of the paper. Needs Python 3.9 or later and nothing else.
# mirror_clause.py: the mirror clause on Hodge data (511 checks).
# yukawa_charges.py: the selection rules on the eleven split data of Section 16 (45 checks);
#   it rewrites yukawa_charges_results.json.
set -euo pipefail
cd "$(dirname "$0")"
python3 -I mirror_clause.py
python3 -I yukawa_charges.py
