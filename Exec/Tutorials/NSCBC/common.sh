# Shared helpers for the NSCBC tutorial rungs.  Source it from a rung's
# run.sh; nothing here is rung-specific.
#
# Environment knobs (all optional):
#   PELE_HOME   PeleC checkout (default: three levels up from this file)
#   USE_MPI     TRUE (default) or FALSE
#   NPROCS      MPI ranks per run (default 4)
#   MAKE_JOBS   parallel compile jobs (default 8)
#   MPIEXEC     launcher (default mpiexec)

TUT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PELE_HOME=${PELE_HOME:-$(cd "$TUT_DIR/../../.." && pwd)}
USE_MPI=${USE_MPI:-TRUE}
NPROCS=${NPROCS:-4}
MAKE_JOBS=${MAKE_JOBS:-8}
MPIEXEC=${MPIEXEC:-mpiexec}

# build_case <Exec/RegTests case> <DIM> [extra make args...]
# Builds the case (and the SUNDIALS TPL once) and sets EXE to the binary.
build_case() {
  local case=$1 dim=$2
  shift 2
  local dir=$PELE_HOME/Exec/RegTests/$case
  mkdir -p runs
  local log=$PWD/runs/build-$case-${dim}d.log
  echo "== building $case (DIM=$dim, USE_MPI=$USE_MPI)"
  (cd "$dir" &&
    make TPL DIM="$dim" USE_MPI="$USE_MPI" "$@" > "$log" 2>&1 &&
    make -j "$MAKE_JOBS" DIM="$dim" USE_MPI="$USE_MPI" "$@" >> "$log" 2>&1) ||
    { echo "build failed: see $log"; exit 1; }
  if [ "$USE_MPI" = "TRUE" ]; then
    EXE=$(ls "$dir"/PeleC${dim}d.*.MPI*.ex 2>/dev/null | head -1)
  else
    EXE=$(ls "$dir"/PeleC${dim}d.*.ex 2>/dev/null | grep -v '\.MPI' | head -1)
  fi
  [ -x "$EXE" ] || { echo "no executable found in $dir"; exit 1; }
  CASE_DIR=$dir
}

# run_variant <name> <inputs file> [ParmParse overrides...]
# Runs EXE in runs/<name> with the given inputs; skips a run that already
# finished (delete its directory to force a rerun).
run_variant() {
  local name=$1 inputs=$2
  shift 2
  local out=runs/$name
  if [ -f "$out/done" ]; then
    echo "   $name: already run (rm -r $out to redo)"
    return
  fi
  mkdir -p "$out"
  echo "   $name"
  (
    cd "$out" || exit 1
    # Success is judged from the log, not the exit code: some compilers
    # (GCC 15) crash in a PelePhysics static destructor after the run has
    # finished and written everything (see the tutorial README).
    if [ "$USE_MPI" = "TRUE" ]; then
      $MPIEXEC -n "$NPROCS" "$EXE" "$inputs" "$@" > run.log 2>&1 || true
    else
      "$EXE" "$inputs" "$@" > run.log 2>&1 || true
    fi
    if grep -q "Ending run" run.log; then
      touch done
    else
      echo "   $name FAILED: see $out/run.log"
    fi
  )
}
