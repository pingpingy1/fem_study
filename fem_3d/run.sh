#!/bin/sh

cd nodes
(set -x; time python3 format_xyz.py)
echo
echo
echo
(set -x; time ./tetgen -cq1.2/10 nodes)
echo
echo
echo
cd ..
(set -x; time gcc -o hamiltonian hamiltonian.c)
echo
echo
echo
(set -x; time ./hamiltonian)
echo
echo
echo
(set -x; time python3 solver.py)

