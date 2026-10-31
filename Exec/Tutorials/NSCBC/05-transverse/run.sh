#!/bin/bash
# Rung 5: oblique waves and vortices -- the transverse weight beta.
# About five minutes on 4 ranks.
set -e
cd "$(dirname "$0")"
source ../common.sh

build_case NSCBC-COVO 2

# (A) a circular pulse in a square with four characteristic faces: the
# front meets the faces at every angle.
PULSE=$CASE_DIR/nscbc-pulse2d.inp
run_variant pulse_hard "$PULSE" pelec.bc_nscbc=0
for b in 1.0 0.8 0.5 0.0; do
  run_variant pulse_b$b "$PULSE" pelec.bc_nscbc_beta=$b
done

# (B) a vortex carried through a characteristic inflow and outflow at
# M = 0.2, and the same vortex in a fully periodic box (no boundary at all:
# the discretisation floor).
COVO=$CASE_DIR/nscbc-covo.inp
run_variant covo_hard "$COVO" pelec.bc_nscbc=0 amr.plot_int=3000
for b in 1.0 0.5 0.2 0.0; do
  run_variant covo_b$b "$COVO" pelec.bc_nscbc_beta=$b amr.plot_int=3000
done
run_variant covo_periodic "$COVO" pelec.bc_nscbc=0 amr.plot_int=3000 \
  geometry.is_periodic="1 1" pelec.lo_bc="Interior Interior" \
  pelec.hi_bc="Interior Interior"
echo "now: python plot.py"
