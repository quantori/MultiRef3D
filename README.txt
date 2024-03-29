# MultiRef3D
This repository hosts files needed to demonstrate and reproduce results presented in the manuscript "A multi-reference poly-conformational method for in silico design, optimization, and repositioning of pharmaceutical compounds illustrated for selected SARS-CoV-2 ligands" by V.Alexandrov, A.Kirpich, O.Kantidze and Y.Gankin

The code was tested under Python 3.7 environment requiring to install the following packages (via conda install):
1. RDkit
2. ODDT

run_score_universe.py scores each compound from a fingerprinted library (either ChEMBL or Enamine) against a set of reference compounds. The output is a list of hits sorted by their total conformer overlap score (Wall) as described in the paper.

proc_fp*.py scripts create fingerprinted universe for screening (either ChEMBL or Enamine).

utils_*.py contain data management and conformer generation / fingerprinting routines.
