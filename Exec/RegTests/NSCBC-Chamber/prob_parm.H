#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real pamb = 1013250.0;

  // Ignition kernel: a half-disc of burnt products centred on the closed-end
  // wall at (problo[0] + kernel_x, kernel_y * Ly), radius kernel_r [cm].  The
  // 1-D flame profile is wrapped around its edge by radial signed distance.
  // kernel_x = 0 is the domain wall (the duct variants); the EB box variants
  // put the closed end at the chamber's EB back wall.
  amrex::Real kernel_r = 0.1;
  amrex::Real kernel_x = 0.0;
  amrex::Real kernel_y = 0.5;

  // Coordinate of the flame in the 1-D PMF profile (see NSCBC-FlameOutflow).
  amrex::Real pmf_flame_loc = 3.0;

  // Laminar flame speed, from the profile's inlet velocity; banner only.
  amrex::Real s_L = 0.0; // filled by amrex_probinit

  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> L = {{1.0}};
  int pmf_N = 0;
  int pmf_M = 0;
  bool pmf_do_average = false;
  amrex::Real* d_pmf_X = nullptr;
  amrex::Real* d_pmf_Y = nullptr;
};

struct ProbParmHost
{
  amrex::Vector<amrex::Real> h_pmf_X;
  amrex::Vector<amrex::Real> h_pmf_Y;
  amrex::Gpu::DeviceVector<amrex::Real> pmf_X;
  amrex::Gpu::DeviceVector<amrex::Real> pmf_Y;
  ProbParmHost() : pmf_X(0), pmf_Y(0) {}
};

#endif
