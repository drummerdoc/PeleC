#!/bin/bash
# Rung 3: the subsonic inflow -- relax_u, forcing, and the NRI register.
# About five minutes on 4 ranks.
set -e
cd "$(dirname "$0")"
source ../common.sh

build_case NSCBC-Acoustic 2
INP=$CASE_DIR/nscbc-acoustic-duct.inp
# The duct: characteristic inflow at x = 0 (target u0), hard ambient
# pressure at x = L = 10 cm (a fully reflecting far end).  t_a = 2 L / c.

# (A) a step in the inflow target: the duct starts at 2000 cm/s, the inlet
# asks for 2100 from t = 0.  20 t_a; the data log gives the mean velocity.
STEP="prob.u_init=2000 prob.u0=2100 prob.force_amp=0 stop_time=0.01152
  amr.plot_files_output=0 amr.data_log=datlog pelec.sum_interval=10"
for k in 0.1 0.3 1 3 10; do
  run_variant step_k$k "$INP" pelec.bc_nscbc_relax_u=$k $STEP
done

# (B) a forced inlet at 0.8 f0 (f0 = 865.1 Hz, the duct's quarter-wave
# frequency): 6 t_a to settle, then 4 forcing periods sampled 24 times each.
FORCE="prob.force_freq=692.1 stop_time=9.24e-3 amr.plot_per=6.02e-5
  amr.plot_vars=density amr.derive_plot_vars=pressure x_velocity"
for k in 0.5 2 5; do
  run_variant force_k${k}_plain "$INP" pelec.bc_nscbc_relax_u=$k $FORCE
  run_variant force_k${k}_ff "$INP" pelec.bc_nscbc_relax_u=$k $FORCE \
    prob.force_feedforward=1
  run_variant force_k${k}_ffnri "$INP" pelec.bc_nscbc_relax_u=$k $FORCE \
    prob.force_feedforward=1 pelec.bc_nscbc_nri=1
done
echo "now: python plot.py"
