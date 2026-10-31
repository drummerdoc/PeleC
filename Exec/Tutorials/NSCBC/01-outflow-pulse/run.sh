#!/bin/bash
# Rung 1: a pressure pulse leaving through a hard-pressure outflow and
# through the characteristic one.  About a minute on 4 ranks.
set -e
cd "$(dirname "$0")"
source ../common.sh

build_case NSCBC-Acoustic 2
INP=$CASE_DIR/nscbc-acoustic.inp
# A plotfile every 20 steps builds the x-t diagram.
OUT="amr.plot_int=20 amr.derive_plot_vars=pressure x_velocity"

run_variant hard "$INP" pelec.bc_nscbc=0 $OUT
run_variant nscbc "$INP" pelec.bc_nscbc=1 pelec.bc_nscbc_sigma=0.25 $OUT
echo "now: python plot.py"
