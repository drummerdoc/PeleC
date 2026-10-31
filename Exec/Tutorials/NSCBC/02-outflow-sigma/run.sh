#!/bin/bash
# Rung 2: the outflow dial sigma -- reflection against holding the mean.
# About five minutes on 4 ranks.
set -e
cd "$(dirname "$0")"
source ../common.sh

build_case NSCBC-Acoustic 2
INP=$CASE_DIR/nscbc-acoustic.inp

# (A) the rung-1 pulse, final state only.
for s in 0 0.25 1 4 16; do
  run_variant pulse_s$s "$INP" pelec.bc_nscbc_sigma=$s amr.plot_int=100000
done
run_variant pulse_s16_ndnr "$INP" pelec.bc_nscbc_sigma=16 \
  pelec.bc_nscbc_ndnr=1 amr.plot_int=100000

# (B) a uniform 0.1% overpressure at rest in the same duct (wall at x = 0,
# outflow at x = L): the "pulse" made wider than the domain.  40 acoustic
# transits; the data log records the domain integrals every 10 steps.
RING="prob.pulse_type=1 prob.width=1.0e5 stop_time=0.01158 max_step=200000
  amr.plot_files_output=0 amr.data_log=datlog pelec.sum_interval=10"
for s in 0 0.25 1 4 16; do
  run_variant ring_s$s "$INP" pelec.bc_nscbc_sigma=$s $RING
done
run_variant ring_s16_ndnr "$INP" pelec.bc_nscbc_sigma=16 \
  pelec.bc_nscbc_ndnr=1 $RING
echo "now: python plot.py"
