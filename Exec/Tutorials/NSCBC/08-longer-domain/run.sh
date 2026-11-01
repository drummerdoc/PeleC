#!/bin/bash
# Rung 8: when the answer is a longer domain -- a vented H2/air chamber with
# its characteristic boundary at the vent, and with the duct extended two
# chamber lengths (air beyond the vent) so the boundary sits in the far
# field; each layout at sigma = 0.25 and 16.  The heaviest rung: about half
# an hour per run on a shared node at 8 (vent) and 24 (plenum) ranks; the
# grids split into 16 / 48 boxes, so give it ranks (NPROCS=16 ./run.sh).
set -e
cd "$(dirname "$0")"
source ../common.sh

build_case NSCBC-Chamber 2
C="prob.pmf_datafile=$CASE_DIR/LiDryer_H2_p1_phi0_4000tu0300.dat
  amr.max_grid_size=12 amr.derive_plot_vars=pressure Temp x_velocity"

# The boundary AT the vent (x = 1.2 cm), with the flame-crossing settings,
# at the default sigma and at sigma = 16.
run_variant vent "$CASE_DIR/chamber-vent.inp" $C
run_variant vent_s16 "$CASE_DIR/chamber-vent.inp" $C pelec.bc_nscbc_sigma=16
# The same chamber with the duct continued to x = 3.6 cm, holding air beyond
# the vent (prob.charge_end): identical boundary settings, two chamber lengths
# further away -- and at sigma = 16, to see whether the far boundary's
# setting still matters.
run_variant plenum "$CASE_DIR/chamber-plenum.inp" $C prob.charge_end=1.2
run_variant plenum_s16 "$CASE_DIR/chamber-plenum.inp" $C prob.charge_end=1.2 \
  pelec.bc_nscbc_sigma=16
echo "now: python3 plot.py"
