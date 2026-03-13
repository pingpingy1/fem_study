#!/bin/sh

python3 create_nodes.py
./triangulate
./hamiltonian
python3 solver.py
