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
  amrex::Real x0 = 0.35;         // pulse centre, fraction of the domain
  amrex::Real width = 0.04;      // pulse width, fraction of the domain
};

struct ProbParmHost
{
  ProbParmHost() = default;
};
#endif
