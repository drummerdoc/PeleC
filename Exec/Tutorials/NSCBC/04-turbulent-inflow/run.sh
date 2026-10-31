#!/bin/bash
# Rung 4: a turbulent inflow through the characteristic inlet, 3-D.
# About five minutes on 4 ranks, plus the 3-D build.
set -e
cd "$(dirname "$0")"
source ../common.sh

# The turbulence: a synthetic HIT box (make_hit.py, numpy only), converted
# to a TurbInflow file by PelePhysics' generator.
GEN_DIR=$PELE_HOME/Submodules/PelePhysics/Support/TurbInflowGenerator
mkdir -p runs
if [ ! -d runs/HIT ]; then
  echo "== building the TurbInflow generator"
  (cd "$GEN_DIR" && make -j "$MAKE_JOBS" DIM=3 USE_MPI=FALSE > /dev/null 2>&1)
  GEN=$(ls "$GEN_DIR"/PeleTurb3d.*.ex | head -1)
  python3 make_hit.py 32 2 runs/hit.dat
  (cd runs && "$GEN" type=turb_box hit_file=hit.dat input_ncell=32 \
    ofile=HIT > hit.log 2>&1)
fi

build_case NSCBC-Acoustic 3
INP=$PWD/duct3d.inp
run_variant defaults "$INP" pelec.bc_nscbc_relax_u=2 pelec.bc_nscbc_relax_t=0.2
run_variant u20 "$INP" pelec.bc_nscbc_relax_u=20 pelec.bc_nscbc_relax_t=0.2
run_variant t20 "$INP" pelec.bc_nscbc_relax_u=2 pelec.bc_nscbc_relax_t=20
run_variant ut20 "$INP" pelec.bc_nscbc_relax_u=20 pelec.bc_nscbc_relax_t=20
run_variant ut20_nri "$INP" pelec.bc_nscbc_relax_u=20 \
  pelec.bc_nscbc_relax_t=20 pelec.bc_nscbc_nri=1
echo "now: python plot.py"
