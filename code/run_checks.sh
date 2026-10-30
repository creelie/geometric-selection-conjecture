#!/usr/bin/env bash
# Runs the checks of the paper. Needs Python 3.9 or later; extension_branch.py also needs sympy.
# mirror_clause.py: the mirror clause on Hodge data (511 checks).
# yukawa_charges.py: the selection rules on the eleven split data of Section 16 (45 checks);
#   it rewrites yukawa_charges_results.json.
# extension_branch.py: the extension branch of the eleven split data, Section 16.3 (61 checks);
#   it rewrites extension_branch_results.json.
# gamma_branch.py: the eleven data and the extension branch on the quotient X_pi / Gamma,
#   Sections 16.2 and 16.4 (152 checks); it rewrites gamma_branch_results.json.
# filtered_branch.py: bundles filtered by the line bundles of the five Gamma-data, Section 16.5
#   (644 checks); it rewrites filtered_branch_results.json.
set -euo pipefail
cd "$(dirname "$0")"
python3 -I mirror_clause.py
python3 -I yukawa_charges.py
python3 -I extension_branch.py
python3 -I gamma_branch.py
python3 -I filtered_branch.py
