#!/bin/sh

PYTHONFILE="format_xyz.py"

while getopts "hc" OPTION; do
    case $OPTION in
        h)
            echo "Usage:"
            echo "$0 [-hc]"
            echo ""
            echo "    -h  display this message"
            echo "    -c  use the cube example"
            echo ""
            exit 0
            ;;
        c)
            PYTHONFILE="create_nodes.py"
            ;;
    esac
done

cd nodes
(set -x; time python3 $PYTHONFILE)
echo
echo
echo
(set -x; time ./tetgen -c nodes)
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

