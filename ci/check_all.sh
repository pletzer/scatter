#!/bin/bash
# Build the compiled variants and check that every variant reproduces the
# reference checksum on a small problem. Run from anywhere; uses the active python.
set -euo pipefail

# usage: check_all.sh [dir]   dir = exercises (default) or solutions
REPO="$(cd "$(dirname "$0")/.." && pwd)"
ROOT="$REPO/${1:-.}"
ROOT="$(cd "$ROOT" && pwd)"
ARGS="-checksum -nx 16 -ny 16 -nc 32"
REF=85.2309      # reference |scat|^2 sum for the small problem
RTOL=1e-4
MPIRUN="${MPIRUN:-mpirun --oversubscribe -n 2}"

check() {  # check <name> <output>
    python - "$1" "$2" "$REF" "$RTOL" <<'PY'
import sys
name, out, ref, rtol = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
val = float(out.strip().splitlines()[-1].split(':')[-1])
ok = abs(val - ref) <= rtol * ref
print(f"{name:10s} {val:.6f}  {'OK' if ok else 'FAIL (expected ' + str(ref) + ')'}")
sys.exit(0 if ok else 1)
PY
}

for d in cext numba openmp; do
    (cd "$ROOT/$d" && rm -rf build && python setup.py build > /dev/null)
done

for d in original vect multiproc numba cext openmp; do
    out=$(cd "$ROOT/$d" && python scatter.py $ARGS)
    check "$d" "$out"
done

# the jax example is complete, so it has no separate solution
if [ -d "$ROOT/jax" ]; then
    out=$(cd "$ROOT/jax" && python scatter.py $ARGS)
    check jax "$out"
fi

out=$(cd "$ROOT/mpi" && $MPIRUN python scatter.py $ARGS)
check mpi "$out"

# -save writes the VTK files (the mpi exercise is unfinished, so only check
# it for the solutions)
saves="vect"
[ "$(basename "$ROOT")" = solutions ] && saves="vect mpi"
for d in $saves; do
    tmp=$(mktemp -d)
    if [ "$d" = mpi ]; then run="$MPIRUN python"; else run="python"; fi
    (cd "$ROOT/$d" && $run scatter.py -nx 8 -ny 8 -nc 16 -save > /dev/null && mv scatter_*.vtk "$tmp")
    test "$(ls "$tmp"/scatter_*.vtk | wc -l)" -eq 20 && echo "save $d   OK"
    rm -rf "$tmp"
done
