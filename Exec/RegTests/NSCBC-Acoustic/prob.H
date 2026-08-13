#ifndef PROB_H
#define PROB_H

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemSpecificFunctions.H"
#include "prob_parm.H"

// ---------------------------------------------------------------------------
//  NSCBC-Acoustic
//
//  A right-running isentropic acoustic pulse in a quiescent gas, launched
//  toward a non-reflecting subsonic outflow at x-hi.  This is the AMReX-side
//  counterpart of check C4 in Verification/NSCBC1D: the pulse should leave the
//  domain and the residual left behind should be a small fraction of the
//  incident amplitude.
//
//  Run it with pelec.bc_nscbc = 0 to see what a hard boundary does instead.
// ---------------------------------------------------------------------------
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x =
    prob_lo[0] + (static_cast<amrex::Real>(i) + 0.5) * dx[0];
  const amrex::Real L = prob_hi[0] - prob_lo[0];

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[0] = 1.0;

  const amrex::Real xc = prob_lo[0] + prob_parm.x0 * L;
  const amrex::Real w = prob_parm.width * L;
  const amrex::Real dp =
    prob_parm.amp * prob_parm.p_amb * std::exp(-std::pow((x - xc) / w, 2));

  // Isentropic, right-running: dp = rho c du, drho = dp / c^2
  const amrex::Real p = prob_parm.p_amb + dp;
  const amrex::Real rho =
    prob_parm.rho_amb + dp / (prob_parm.cs_amb * prob_parm.cs_amb);
  const amrex::Real u = dp / (prob_parm.rho_amb * prob_parm.cs_amb);

  amrex::Real T = 0.0, eint = 0.0;
  eos.RYP2T(rho, massfrac, p, T);
  eos.RTY2E(rho, T, massfrac, eint);

  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = rho * eint;
  state(i, j, k, UEDEN) = rho * (eint + 0.5 * u * u);
  state(i, j, k, UTEMP) = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) = rho * massfrac[n];
  }
}

// The ordinary ghost fill.  Only reached when pelec.bc_nscbc = 0, or for
// points where bcnormal_nscbc returns `off`.  Deliberately crude -- a hard
// imposition of the ambient state -- so that the comparison against the
// characteristic treatment is stark.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real* s_int,
  amrex::Real* s_ext,
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[0] = 1.0;

  const amrex::Real rho = prob_parm.rho_amb;
  const amrex::Real u = s_int[UMX] / s_int[URHO];
  amrex::Real T = 0.0, eint = 0.0;
  eos.RYP2T(rho, massfrac, prob_parm.p_amb, T);
  eos.RTY2E(rho, T, massfrac, eint);

  s_ext[URHO] = rho;
  s_ext[UMX] = rho * u;
  s_ext[UMY] = 0.0;
  s_ext[UMZ] = 0.0;
  s_ext[UEINT] = rho * eint;
  s_ext[UEDEN] = rho * (eint + 0.5 * u * u);
  s_ext[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    s_ext[UFS + n] = rho * massfrac[n];
  }
}

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
  // The characteristic treatment.  For a pure non-reflecting outflow the
  // target pressure is the only thing that needs to be specified.
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static pc_nscbc::Target bcnormal_nscbc(
    const amrex::Real* /*x[AMREX_SPACEDIM]*/,
    const amrex::Real* /*s_int[NVAR]*/,
    const int idir,
    const int sgn,
    const amrex::Real /*time*/,
    amrex::GeometryData const& /*geomdata*/,
    ProbParmDevice const& prob_parm)
  {
    pc_nscbc::Target t;
    if ((idir == 0) && (sgn == -1)) {
      t.type = pc_nscbc::Type::outflow;
      t.p = prob_parm.p_amb;
    }
    return t;
  }
};
using ProblemSpecificFunctions = MyProblemSpecificFunctions;

void pc_prob_close();

#endif
