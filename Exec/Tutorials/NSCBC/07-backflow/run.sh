#!/bin/bash
# Rung 7: flow entering through an outflow -- reversal, and sustained backflow.
# About three minutes on 4 ranks.
set -e
cd "$(dirname "$0")"
source ../common.sh

build_case NSCBC-Acoustic 2
INP=$CASE_DIR/nscbc-acoustic.inp
# The flush test: a 10 cm duct, both ends characteristic outflows with an
# ambient (300 K air) reservoir behind them, filled with gas at 600 K and at
# rest.  The x = 0 target is 0.5% below ambient -- a sink -- so a stream is
# drawn in through the x = L outflow and out through x = 0.
FLUSH="prob.pulse_type=1 prob.amp=0.0 prob.T_init=600.0
  prob.p_lo_offset=-5066.0 pelec.lo_bc=UserBC Interior amr.n_cell=200 8
  stop_time=0.03 max_step=1000000 amr.plot_int=-1 amr.plot_per=5.0e-3
  amr.data_log=datlog pelec.sum_interval=20"
run_variant frozen "$INP" $FLUSH pelec.bc_nscbc_backflow_material=0
run_variant reservoir "$INP" $FLUSH pelec.bc_nscbc_backflow_material=1
echo "now: python3 plot.py"
