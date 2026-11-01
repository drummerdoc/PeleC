#!/bin/bash
# Rung 6: what may cross an outflow -- a hot front with and without heat
# release leaving through the boundary.  The heaviest rung: each flame run
# is about 35 min on 4 ranks; the case splits into 12 boxes, so up to 12
# ranks help (NPROCS=12 ./run.sh).
set -e
cd "$(dirname "$0")"
source ../common.sh

build_case NSCBC-FlameOutflow 2
INP=$CASE_DIR/nscbc-flameexit.inp

# The exit test: a wrinkled H2/air flame (phi = 0.4) carried out through the
# outflow at 40 S_L, inflow hard (prob.nscbc_inflow = 0 in the inputs).
# extrap_material stays off: it is for fronts that sit, not ones that leave.
C="prob.pmf_datafile=$CASE_DIR/LiDryer_H2_p1_phi0_4000tu0300.dat
  pelec.bc_nscbc_extrap_material=0 amr.max_grid_size=16 amr.plot_per=1.0e-5
  amr.derive_plot_vars=pressure x_velocity"
CLOSURES="pelec.bc_nscbc_extrap_temperature=1 pelec.bc_nscbc_beta_s=0"
DEFAULTS="pelec.bc_nscbc_extrap_temperature=0 pelec.bc_nscbc_beta_s=1"

# (A) An inert hot front: the same temperature and composition profile with
# the chemistry off, so only entropy and composition cross.
run_variant inert_hard "$INP" $C pelec.bc_nscbc=0 pelec.do_react=0
run_variant inert_defaults "$INP" $C $DEFAULTS pelec.do_react=0
run_variant inert_closures "$INP" $C $CLOSURES pelec.do_react=0
run_variant inert_s16_defaults "$INP" $C $DEFAULTS pelec.do_react=0 \
  pelec.bc_nscbc_sigma=16
run_variant inert_s16_closures "$INP" $C $CLOSURES pelec.do_react=0 \
  pelec.bc_nscbc_sigma=16
# (B) The flame.
run_variant flame_hard "$INP" $C pelec.bc_nscbc=0
run_variant flame_defaults "$INP" $C $DEFAULTS
run_variant flame_closures "$INP" $C $CLOSURES
run_variant flame_s16_defaults "$INP" $C $DEFAULTS pelec.bc_nscbc_sigma=16
run_variant flame_s16_closures "$INP" $C $CLOSURES pelec.bc_nscbc_sigma=16
run_variant flame_s16_closures_ndnr "$INP" $C $CLOSURES \
  pelec.bc_nscbc_sigma=16 pelec.bc_nscbc_ndnr=1
echo "now: python3 plot.py"
