#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p_amb = 1.01325e6; // ambient / outflow target [dyn/cm^2]
  amrex::Real T_amb = 300.0;     // ambient temperature [K]
  amrex::Real rho_amb = 0.0;     // filled by amrex_probinit
  amrex::Real cs_amb = 0.0;      // filled by amrex_probinit
  amrex::Real amp = 1.0e-3;      // pulse amplitude, relative to p_amb

  // 0 = PLANAR pulse running along +x into a single characteristic face.
  //     The historical configuration; keeps its regression gold.
  // 1 = RADIAL pulse at the box centre, every face characteristic.  The pulse
  //     reaches the faces at normal incidence, the edges at 45 degrees and the
  //     corners at the body diagonal, all in one run, which is what makes it
  //     the test of corner and edge ownership rather than of the 1-D algebra.
  int pulse_type = 0;

  // Pulse centre and width, as fractions of the domain extent in each
  // direction.  For pulse_type = 0 only the first component of x0 is used and
  // the width is measured against the first direction; for pulse_type = 1 the
  // centre is a point and the width is a single physical length, so that the
  // pulse is a sphere and not an ellipsoid on an anisotropic mesh.
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> x0 = {
    {AMREX_D_DECL(0.35, 0.5, 0.5)}};
  amrex::Real width = 0.04;

  // Filled by amrex_probinit: the absolute pulse centre and the absolute
  // width, so that pc_initdata does no geometry arithmetic and cannot
  // disagree with the banner.
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> xc = {
    {AMREX_D_DECL(0.0, 0.0, 0.0)}};
  amrex::Real w = 1.0;
};

struct ProbParmHost
{
  ProbParmHost() = default;
};
#endif
