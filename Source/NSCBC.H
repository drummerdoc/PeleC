#ifndef NSCBC_H
#define NSCBC_H
// ============================================================================
//  NSCBC.H -- characteristic (Navier-Stokes Characteristic Boundary Condition)
//             ghost-cell fill for subsonic inflow and outflow faces.
//
//  This header is deliberately free of any dependency on PeleC.H so that it can
//  also be compiled by the standalone 1-D verification driver in
//  Verification/NSCBC1D.  Its only dependencies are IndexDefines.H, the
//  PelePhysics EOS, and AMReX.
//
//  ---------------------------------------------------------------------------
//  FORMULATION
//  ---------------------------------------------------------------------------
//  Everything is written in an "outward-normal" frame,
//
//      u_out = n_sgn * u_n ,      n_sgn = -sgn = +1 on a hi face, -1 on a lo
//                                 face,
//
//  so that a single code path serves all six faces and an outflow always has
//  u_out > 0.  In that frame the wave speeds along the outward normal are
//
//      lambda_-  = u_out - c   (always < 0 subsonically -> ALWAYS INCOMING)
//      lambda_0  = u_out       (outgoing at an outflow, incoming at an inflow)
//      lambda_+  = u_out + c   (always > 0 subsonically -> ALWAYS OUTGOING)
//
//  The lambda_0 family carries the entropy, the two tangential velocities, the
//  NUM_SPECIES mass fractions and any advected/auxiliary/passive scalars.  That
//  the lambda_0 family flips direction between inflow and outflow is the whole
//  content of the inflow/outflow distinction:
//
//      subsonic OUTFLOW: 1 incoming characteristic  -> specify 1 quantity (p)
//      subsonic INFLOW : D+NUM_SPECIES+1 incoming   -> specify u, T, Y
//
//  We work with the linearised characteristic *invariants* rather than with the
//  wave amplitudes L_i and their normal derivatives:
//
//      R_+ = u_out + p/(rho c)      (outgoing)
//      R_- = u_out - p/(rho c)      (incoming)
//      S   = rho - p/c^2            (entropy, linearised)
//
//  with the acoustic impedance rho*c FROZEN at the boundary cell N.  Freezing
//  is not an optimisation: R_+ is only an invariant of the system linearised
//  about a fixed state, so if rho*c were re-evaluated per stencil cell then
//  R_+(N) - R_+(N-1) would be a difference of two different variables and its
//  use as a slope would be meaningless.
//
//  Outgoing quantities (R_+ always; S or T, u_t, Y, scalars at an outflow) are
//  extrapolated outward from the interior with a minmod-limited linear slope.
//  Which of S or T is extrapolated is Params::extrap_temperature -- see there;
//  it matters because the diffusion operator reads these same ghost cells.
//  Incoming quantities are set from a target state through the classical
//  Poinsot-Lele relaxation:
//
//      L_incoming = K (phi - phi_target) ,   K = coeff * (1 - M^2) * c / L_ref
//
//  converted to an increment on the invariant over the ghost offset,
//
//      R_-,ghost = R_-,N + layer * dx * dR_-/dn ,
//      dR_-/dn   = -L_incoming / (lambda_- * rho c) .
//
//  K is a RATE [1/s] built from a physical length, so the boundary relaxation
//  time tau = 1/K is independent of the mesh and of the ghost layer count, and
//  literature values of the relaxation coefficient transfer directly.  The
//  boundary pressure then relaxes toward the far field with
//
//      tau_relax = L_ref / (sigma (1 - M^2) c)
//
//  which at sigma = 0.25 and low Mach is about four acoustic transit times.
//  That is the entire content of "sigma = 0.25 is a good default".
//
//  Two notes on what this is NOT:
//
//    1. Transverse terms are available, weighted by Params::beta, and the
//       reaction source is available, weighted by Params::beta_s -- see the
//       guidance on those members.  Both are computed cell-locally from the
//       boundary cell and its tangential neighbours, so neither needs
//       sources_for_hydro (which does not exist during FillPatch).  The
//       VISCOUS source is genuinely absent.  With beta = beta_s = 1 (the
//       defaults) both are discarded and this is the pure 1-D-normal LODI
//       limit.
//
//       Neither correction rescues an outflow placed on a flame.  Measured in
//       Exec/RegTests/NSCBC-FlameOutflow, beta_s moves the mean pressure error
//       there by under 1%, because the dominant error is not the heat release:
//       it is a ghost-pressure bias driven by the normal velocity gradient,
//       dp ~ rho c L_ref (du_out/dn) / sigma, which is present with chemistry
//       switched off entirely.  Place outflows where the flow is smooth.
//
//       That the transverse terms fit here at all is a consequence of writing
//       the relaxation as a GRADIENT increment: they are additive on
//       L_incoming, so there is a slot for them.  A hard far-field pin of R_-
//       has no such slot, which is why pin_farfield ignores beta.
//
//    2. It is not a time-integrated boundary condition.  The fill is a pure
//       algebraic function of the current interior state, so it is idempotent
//       under the repeated FillPatch calls of an SDC iteration.
//
//  ---------------------------------------------------------------------------
//  SIGN CONVENTIONS
//  ---------------------------------------------------------------------------
//  Every user-facing relaxation coefficient (sigma, relax_u, relax_t) is
//  POSITIVE.  The legacy Fortran required relax_T to be negative, and its
//  documentation and its code disagreed about the sign of sigma_out; all such
//  sign handling is internal here.  Params::validate() rejects negatives.
//
//  ---------------------------------------------------------------------------
//  GPU PORTABILITY
//  ---------------------------------------------------------------------------
//  Finiteness guards use amrex::Math::isfinite, NOT std::isfinite. AMReX_Math.H
//  documents why: SYCL's fast-math makes std::isfinite evaluate to false even
//  at -O1.  Every use of it in this file is a guard of the form
//
//      if (!isfinite(x)) { fall back; bump(counter); }
//
//  so under a SYCL build with std::isfinite EVERY guard would fire and the
//  boundary would silently degrade to a zero-gradient copy -- correct-looking
//  code producing a quietly different boundary condition on one vendor's GPUs.
//  The diagnostic counters would at least make it visible, which is the second
//  reason they exist.
//
//  There is no dynamic allocation, no virtual dispatch, no recursion, and no
//  printf or Abort on the device path; the fallback counters are incremented
//  with amrex::HostDevice::Atomic::Add.
//
//  The one cost worth knowing about is stack footprint.  The caller holds
//  4*NVAR doubles of stencil and ghost state plus 2*(D-1)*NVAR of tangential
//  neighbours, and apply() adds three CellPrim plus the species work arrays --
//  order 10*NVAR doubles per thread in 3-D, which is 2.2 kB for a 21-species
//  mechanism and 4.8 kB for a 53-species one.  That exceeds the 255-register
//  budget and will spill to local memory.  It is tolerable because the launch
//  is over a thin ghost region rather than the whole domain, but it is the
//  first thing to look at if a boundary fill ever shows up in a GPU profile.
//
//  ---------------------------------------------------------------------------
//  EOS PORTABILITY
//  ---------------------------------------------------------------------------
//  Only density-carrying ("R...") EOS entry points are used:
//      RTY2P, RTY2Cs, RYP2T, RTY2E, REY2T, PYT2RE
//  The shortcut forms TY2Cp, TY2G, EY2T, RG2dpde, RPE2dpdr_e, T2Hi et al.
//  amrex::Error() out under the SRK real-gas EOS.  Keeping to the list above
//  makes this kernel real-gas capable at no extra cost.  Note also that SRK's
//  RYP2E is implemented as RYP2T followed by RTY2E, so calling those two
//  directly costs one Newton solve instead of two.
// ============================================================================

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Math.H>
#include <cmath>

#include "IndexDefines.H"
#include "Constants.H"
#include "PelePhysics.H"

namespace pc_nscbc {

// ---------------------------------------------------------------------------
//  Boundary type, decided per boundary POINT by the problem's bcnormal_nscbc().
//  A single face may carry a mixture (a jet inlet surrounded by outflow).
// ---------------------------------------------------------------------------
enum struct Type : int { off = 0, outflow = 1, inflow = 2 };

// ---------------------------------------------------------------------------
//  Numerical parameters.  POD, captured host-side into the boundary functor;
//  never read ParmParse or a class static from inside a device kernel.
// ---------------------------------------------------------------------------
struct Params
{
  // Outflow pressure relaxation, classical Poinsot-Lele sigma:
  //     K = sigma * (1 - M^2) * c / L_ref      [1/s]
  // 0 is perfectly non-reflecting and leaves the mean pressure unanchored;
  // the literature band is 0.15-0.3 (Rudy & Strikwerda 1980 give ~0.27 as
  // optimal for a 1-D pulse).
  amrex::Real sigma = 0.25;

  // Inflow normal-velocity relaxation.  Larger imposes the target velocity
  // more strongly and reflects more; >~10 is a hard Dirichlet in disguise.
  amrex::Real relax_u = 2.0;

  // Inflow temperature and tangential-velocity relaxation.
  amrex::Real relax_t = 0.2;

  // Weight on the transverse terms in the modelled incoming acoustic:
  //     L_in = K (phi - phi_target) - (1 - beta) * T_transverse
  //   1  discard them entirely -- the 1-D-normal LODI limit
  //   0  include them fully
  //   <0 use the local Mach number, |M|, clamped to [0,1]
  // beta = 1 reproduces the Phase 1 behaviour bit for bit, which is why it is
  // the default: enabling the transverse terms is an opt-in change.
  //
  // Guidance, measured in Exec/RegTests/NSCBC-COVO rather than assumed.
  //
  // For a plane wave meeting the boundary at angle theta to the normal, the
  // obliqueness error that the transverse term has to cancel works out to
  //
  //     (1 - beta)_opt = cos(theta) / (1 + cos(theta))
  //     =>  beta_opt   = 1 - cos(theta)/(1 + cos(theta))
  //
  // which is 0.5 at normal incidence and rises toward 1 as the wave becomes
  // grazing (0.67 at 60 deg, 0.79 at 75 deg).  Both measured optima land on
  // that curve: the convected vortex, a broad near-normal structure, is best
  // at beta = 0.5; the circular pulse, whose surviving arc crossed the faces
  // near-tangentially, is best at beta = 0.8.
  //
  //     COVO   rms residual vs the no-boundary floor:
  //            hard 14.0x,  beta=1 10.0x,  beta=0.5 3.8x,  beta=0.2 16.1x,
  //            beta=0 132x
  //     pulse  front amplitude spread around the arc:
  //            hard 21.4%,  beta=1 0.90%,  beta=0.8 0.066%,  beta=0 2.9%
  //
  // Note the shape of that curve: too little correction costs a factor of a
  // few, too much is catastrophic.  beta below ~0.3 is worse than no
  // transverse terms at all, and beta = 0 is close to unstable.  This is why
  // the default is 1 (the Phase 1 behaviour) and not the optimum -- a wrong
  // beta is far more dangerous than an absent one.  Start at 0.5, and raise
  // it toward 1 if the boundary is mostly met at grazing incidence.
  //
  // beta < 0 selects the local Mach number, which is the legacy Fortran's
  // convention.  The measurements above do not support it: at M = 0.2 it
  // gives beta = 0.2, which is on the wrong side of the minimum.
  amrex::Real beta = 1.0;

  // Weight on the REACTION source term in the modelled incoming acoustic,
  // with the same convention as beta:
  //     L_in = K (phi - phi_target) - (1 - beta) T - (1 - beta_s) S_reaction
  //   1  discard it -- the default, and the only correct choice when the
  //      boundary is far enough from the flame that heat release there is
  //      negligible
  //   0  include it fully
  //
  // What this corrects.  A flame near an outflow raises the pressure in the
  // boundary cell at a rate the 1-D relaxation has no model for, so sigma has
  // to absorb it and the mean pressure sits off target.  S_reaction supplies
  // the missing dp/dt.
  //
  // Note that this is NOT obtainable from PeleC's srcq(QPRES): that expression
  // (Hydro.H) is built only from the density and internal-energy sources, and
  // chemistry supplies neither -- mass is conserved and PeleC's constant-volume
  // reactor conserves total energy because the formation enthalpies live inside
  // e.  The entire effect is compositional, dp/dY|_{rho,e}, which srcq(QPRES)
  // does not carry.  See reaction_dpdt() below.
  amrex::Real beta_s = 1.0;

  // Close the lambda_0 family on TEMPERATURE rather than on the linearised
  // entropy invariant S = rho - p/c^2.
  //
  // Both are legitimate: the lambda_0 family carries entropy, and T is a
  // function of (rho, p, Y), so extrapolating either one and deriving the other
  // closes the family.  For a smooth flow they agree to first order.  They do
  // NOT agree across a steep temperature ramp, and the difference matters for a
  // reason that has nothing to do with the hyperbolic problem: PeleC's
  // diffusion operator forms the conductive and species fluxes at a physical
  // boundary face from these ghost cells, so whatever normal temperature
  // gradient the ghost carries IS the heat flux leaving the domain.
  //
  // With the entropy closure, T is whatever the EOS returns from the
  // extrapolated rho and the acoustically-set p, and nothing constrains its
  // gradient.  Check C8 in Verification/NSCBC1D measures the consequence on a
  // flame-like ramp (2e4 K/cm): the ghost overstates the face temperature
  // gradient by 40%.  With this set, T is extrapolated with the same
  // minmod-limited slope as everything else and rho follows from the EOS, so
  // the conductive flux at the face is a controlled continuation of the
  // interior one.
  //
  // Default false: the entropy closure is the historical behaviour and the
  // textbook one, and switching it changes every outflow result, so it is
  // opt-in.  Turn it on when a thermal or compositional structure is anywhere
  // near the boundary.
  bool extrap_temperature = false;

  // Continue the MATERIAL part of the incoming invariant R_- into the ghost
  // at an outflow, so that a dilatational velocity gradient crossing the face
  // is not converted into ghost pressure.
  //
  // The problem this solves.  Linearised about the boundary cell,
  //
  //     outgoing acoustic wave:            dR_+ = 2 du,   dR_- = 0
  //     material structure (dp ~ 0):       dR_+ = dR_- = du
  //
  // The default treatment extrapolates R_+ with its full interior slope and
  // gives R_- only the relaxation increment, so across material structure the
  // ghost picks up p_g - p_N = 1/2 rho c l du_out per layer -- the bias of
  // check C9(a) in Verification/NSCBC1D and the dominant error at a
  // front-crossing outflow (Exec/RegTests/NSCBC-FlameOutflow).  Continuing
  // the material slope of R_- as well makes the ghost carry the interior's
  // u and p slopes exactly and the bias vanishes, at ANY sigma.
  //
  // Why the slope cannot simply be measured from R_-.  The interior gradient
  // of R_- also contains every piece of INCOMING acoustics near the boundary
  // -- including the anchoring wave the relaxation itself launches -- and
  // extrapolating that back into the ghost is a positive feedback loop: it
  // multiplies the effective relaxation rate (measured 16x in the C5 rate
  // check) and reflects an outgoing pulse at 13% instead of 0.8%.  The
  // 1-D acoustic pair alone cannot distinguish material structure from
  // counter-running waves; what does distinguish them is the entropy family,
  // which is OUTGOING at an outflow and carries no acoustic content:
  //
  //     dS = drho - dp/c^2      -> 0        for any acoustic wave
  //                             -> drho     for material structure
  //
  // Steady continuity through the face (d(rho u) = 0, with the momentum-flux
  // pressure gradient dp = -rho u du) turns that into a bound on the material
  // velocity slope,
  //
  //     du_mat = -u_out dS / (rho (1 - M^2)) ,   dR_-^mat = du_mat (1 + M) ,
  //
  // and the applied slope is minmod(measured dR_-, dR_-^mat): across a flame
  // the two agree and the measured slope is used exactly; for pure acoustics
  // the entropy bound is ~0 and the continuation shuts itself off, so the
  // feedback path is closed -- every quantity the incoming model consumes
  // rides an outgoing family.  A structure that advects without dilatation
  // (a contact: du = 0, drho != 0) has measured dR_- = 0 and is equally
  // unaffected.
  //
  // The bound assumes the structure is quasi-steady in the LAB frame; a front
  // advancing against the flow at speed V has its true slope underestimated
  // by u/(u - V), so the continuation is partial there -- conservative, never
  // amplifying.  Keep sigma > 0 with this on: the relaxation increment is
  // still the only anchoring.
  //
  // Measured (Exec/RegTests/NSCBC-FlameOutflow, quasi-frozen front, sigma=1):
  // 5% alone, 42% together with extrap_temperature -- the two closures
  // compound into a fully consistent material continuation.  Measured on the
  // exit test of the same case (a front actually LEAVING at U = 40 S_L):
  // at sigma = 0.25 the continuation turns the crossing into a runaway that
  // drains the domain and quenches the flame -- the quasi-steady model does
  // not describe a fast transit, and the weak relaxation cannot answer for
  // it.  At sigma = 16 it is safe but unhelpful.
  //
  // Default false: it changes outflow results near structures, so it is
  // opt-in, like extrap_temperature.  Turn it on for material structure that
  // SITS near the outflow (a held or slowly drifting flame, a density
  // gradient, a wake); leave it off, or pair it with a large sigma, when a
  // front actually passes out through the face.  pin_farfield ignores it: a
  // value pin has no slot for a slope.
  bool extrap_material = false;

  // Boundary-normal reference length, probhi[idir] - problo[idir].  Filled per
  // face by the caller.  Only the ratio sigma/L_ref is physical, so this is
  // not a user knob; it exists so that sigma keeps its literature meaning.
  amrex::Real L_ref = 1.0;

  // Order of the outgoing-invariant extrapolation:
  //   1 = zeroth-order copy (constant across ghost layers)
  //   2 = minmod-limited linear in the ghost layer index (default)
  // Exposed so that a bit-comparison against the standalone 1-D driver can be
  // done without recompiling.  Verification knob, not a physics knob.
  int order = 2;

  // Replace the Poinsot-Lele relaxation of the incoming acoustic by a hard pin
  // of R_- to the far-field value (u_out = 0, p = p_target).  It is
  // simultaneously non-reflecting and anchored, but it anchors to
  // p_target + rho c u_out rather than to p_target, and it is a value
  // constraint rather than a rate, so its effective relaxation rate is c/dx
  // -- i.e. it does not converge under grid refinement.
  // Right for an open boundary onto a large quiescent reservoir; wrong for a
  // duct exhausting into a plenum whose true mean pressure is not p_target.
  // sigma is ignored when this is set.
  bool pin_farfield = false;

  AMREX_GPU_HOST bool validate(std::string& why) const
  {
    if (sigma < 0.0) {
      why = "bc_nscbc_sigma must be >= 0 (all NSCBC relaxation coefficients "
            "are positive in PeleC; internal signs are handled by the kernel)";
      return false;
    }
    if (relax_u < 0.0) {
      why = "bc_nscbc_relax_u must be >= 0";
      return false;
    }
    if (relax_t < 0.0) {
      why = "bc_nscbc_relax_t must be >= 0";
      return false;
    }
    if (beta_s > 1.0) {
      why = "bc_nscbc_beta_s must be <= 1";
      return false;
    }
    if (beta > 1.0) {
      why = "bc_nscbc_beta must be <= 1 (1 discards the transverse terms, 0 "
            "includes them fully, negative selects the local Mach number)";
      return false;
    }
    if ((order != 1) && (order != 2)) {
      why = "bc_nscbc_order must be 1 or 2";
      return false;
    }
    if (L_ref <= 0.0) {
      why = "NSCBC reference length must be > 0";
      return false;
    }
    return true;
  }
};

// ---------------------------------------------------------------------------
//  Target state returned per boundary point by the problem hook.
//
//  Only the entries relevant to `type` are read:
//    outflow -> p
//    inflow  -> u[], T, Y[]
//  There is deliberately no relaxation coefficient for composition: at a
//  subsonic inflow every species characteristic is incoming, so hard
//  imposition is well posed; at a subsonic outflow every species
//  characteristic is outgoing, so composition is extrapolated and imposing it
//  would over-specify the problem.
// ---------------------------------------------------------------------------
struct Target
{
  Type type = Type::off;
  amrex::Real p = 0.0;                // outflow: target pressure [dyn/cm^2]
  amrex::Real u[3] = {0.0, 0.0, 0.0}; // inflow: target velocity [cm/s]
  amrex::Real T = 0.0;                // inflow: target temperature [K]
  amrex::Real Y[NUM_SPECIES] = {0.0}; // inflow: target mass fractions
};

// ---------------------------------------------------------------------------
//  Tangential neighbours of the BOUNDARY CELL, for the transverse terms.
//
//  The transverse contribution to the incoming acoustic amplitude is, in the
//  outward-normal frame and summed over the tangential directions t,
//
//      T_in = sum_t [ u_t ( dp/dt - rho c du_out/dt ) ]  +  gamma p sum_t
//      du_t/dt
//
//  which reduces to Motheau's T1 on a hi face and to his T4 on a lo face -- the
//  outward frame absorbs the sign flip that makes those two look like different
//  expressions.  Every derivative is evaluated at the boundary cell.
//
//  The caller supplies the two tangential neighbours per direction with their
//  indices ALREADY CLAMPED into the domain, so only interior cells are ever
//  read and the fill stays a pure function of valid data.  At a corner the
//  clamp collapses the stencil, which is why the spacing is passed as an
//  explicit reciprocal rather than assumed: inv_dt = 1/(x_p - x_m), and zero
//  when both neighbours land on the same cell.  That degrades a corner to a
//  one-sided or vanishing transverse derivative instead of reading a ghost
//  cell another thread is writing.
// ---------------------------------------------------------------------------
struct Transverse
{
  const amrex::Real* sm[AMREX_SPACEDIM] = {};
  const amrex::Real* sp[AMREX_SPACEDIM] = {};
  amrex::Real inv_dt[AMREX_SPACEDIM] = {};
  bool valid = false;
};

// ---------------------------------------------------------------------------
//  Diagnostic counters.  Every fallback in this kernel is counted; a silent
//  fallback is a bug that will not be found.  Indices into a small device
//  array reduced at the end of the run.  The storage is amrex::Long, not int:
//  the counters are only reset when they are REPORTED, and with the default
//  sum_interval = -1 they never are -- the supersonic counter on a large 3-D
//  run increments once per supersonic ghost cell per fill and walks past 2^31
//  in a few times 10^4 steps, and an overflowing signed atomic is exactly the
//  kind of quiet corruption these counters exist to prevent.  amrex::Gpu::
//  Atomic::Add has a native Long specialisation on every backend.
// ---------------------------------------------------------------------------
namespace Diag {
enum : int {
  supersonic = 0,  // |M| >= 1: zero-gradient copy is exact, not a failure
  reversed,        // flow ran the wrong way through the face
  body_state,      // an EB-covered cell was found in the stencil
  eos_failure,     // an EOS call returned a non-finite or non-positive value
  floored,         // a positivity floor was hit
  transverse_drop, // a transverse derivative was unavailable or unusable
  source_drop,     // a modelled source term was unavailable or unusable
  structure,       // material structure sat in an outflow boundary cell:
                   // |dS| > 5% of rho per cell.  ADVISORY, not a fallback --
                   // if a front is crossing this outflow, sigma of O(10) is
                   // what survives it, and the sigma = 0.25 default can push
                   // the front back and abort the run (see the flame-exit
                   // table in NSCBC-FlameOutflow/README.md).
  count
};
}

// ===========================================================================
//  Small helpers
// ===========================================================================

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
finite_or(const amrex::Real x, const amrex::Real fallback) noexcept
{
  return amrex::Math::isfinite(x) ? x : fallback;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
minmod(const amrex::Real a, const amrex::Real b) noexcept
{
  if (a * b <= 0.0) {
    return 0.0;
  }
  return (std::abs(a) < std::abs(b)) ? a : b;
}

// True if a state looks like EB body state rather than fluid.  PeleC fills
// covered cells with -1 per component (IO.cpp) or with 0 when
// eb_zero_body_state is set; either way the density is non-positive.  The
// boundary functor has no EBCellFlag, so this predicate is how the stencil
// notices.  Cheap, and it also catches genuinely broken states.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
is_body_state(const amrex::Real s[NVAR]) noexcept
{
  return !(s[URHO] > constants::very_small_num()) ||
         !amrex::Math::isfinite(s[URHO]) || !amrex::Math::isfinite(s[UEDEN]);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
copy_state(const amrex::Real src[NVAR], amrex::Real dst[NVAR]) noexcept
{
  for (int n = 0; n < NVAR; n++) {
    dst[n] = src[n];
  }
}

// ---------------------------------------------------------------------------
//  Primitive decomposition of one stencil cell.
//
//  Temperature is taken from the state rather than recovered with REY2T.  T is
//  a maintained component of PeleC's conserved state (computeTemp runs after
//  every update) and using it keeps every EOS call on this path non-iterative,
//  which matters under SRK.  This is the same choice the existing
//  characteristic-outflow bcnormal implementations make (see EB-C10/prob.H).
//  If T is not usable we fall back to the iterative recovery.
// ---------------------------------------------------------------------------
struct CellPrim
{
  amrex::Real rho, u[3], T, p, c, e;
  amrex::Real Y[NUM_SPECIES];
  bool ok;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE CellPrim
cell_primitives(const amrex::Real s[NVAR]) noexcept
{
  auto eos = pele::physics::PhysicsType::eos();
  CellPrim q;
  q.ok = true;

  q.rho = s[URHO];
  const amrex::Real rhoinv = 1.0 / q.rho;
  q.u[0] = s[UMX] * rhoinv;
  q.u[1] = s[UMY] * rhoinv;
  q.u[2] = s[UMZ] * rhoinv;

  amrex::Real Ysum = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    q.Y[n] = amrex::max<amrex::Real>(s[UFS + n] * rhoinv, 0.0);
    Ysum += q.Y[n];
  }
  if (Ysum > constants::very_small_num()) {
    const amrex::Real inv = 1.0 / Ysum;
    for (int n = 0; n < NUM_SPECIES; n++) {
      q.Y[n] *= inv;
    }
  } else {
    q.ok = false;
  }

  const amrex::Real ke =
    0.5 * (q.u[0] * q.u[0] + q.u[1] * q.u[1] + q.u[2] * q.u[2]);
  q.e = s[UEDEN] * rhoinv - ke;

  q.T = s[UTEMP];
  if (!(q.T > 0.0) || !amrex::Math::isfinite(q.T)) {
    q.T = 300.0; // initial guess for the iterative fallback
    eos.REY2T(q.rho, q.e, q.Y, q.T);
  }

  eos.RTY2P(q.rho, q.T, q.Y, q.p);
  eos.RTY2Cs(q.rho, q.T, q.Y, q.c);

  if (
    !amrex::Math::isfinite(q.p) || !(q.p > 0.0) ||
    !amrex::Math::isfinite(q.c) || !(q.c > 0.0) ||
    !amrex::Math::isfinite(q.T) || !(q.T > 0.0)) {
    q.ok = false;
  }
  return q;
}

// ---------------------------------------------------------------------------
//  reaction_dpdt -- the rate of change of pressure at a cell due to chemistry
//  alone, holding density and total internal energy fixed.
//
//  Chemistry conserves mass and, in PeleC's constant-volume reactor, total
//  internal energy: the formation enthalpies are carried inside e.  So it
//  changes the pressure entirely through the COMPOSITION,
//
//      dp/dt|_react = sum_k (dp/dY_k)|_{rho,e} * wdot_k / rho
//
//  For an ideal-gas mixture that derivative is available in closed form,
//
//      (dp/dY_k)|_{rho,e} = rho Ru T / W_k  -  (p/T)(e_k / c_v)
//
//  giving a mole-change term and the familiar heat-release term:
//
//      dp/dt = Ru T sum_k wdot_k/W_k  -  p/(rho T c_v) sum_k e_k wdot_k
//
//  For a real-gas EOS there is no such one-liner, and SRK does not implement
//  the composition derivative at all.  The fallback exploits the fact that
//  sum_k wdot_k = 0, so the composition perturbation along the reaction path
//  already preserves sum Y = 1: a single DIRECTIONAL finite difference along
//  that path gives exactly the required derivative, at the cost of one extra
//  (rho, e, Y) -> T solve rather than NUM_SPECIES of them.
// ---------------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
reaction_dpdt(CellPrim const& q, bool& ok) noexcept
{
  auto eos = pele::physics::PhysicsType::eos();
  ok = true;

  amrex::Real wdot[NUM_SPECIES];
  eos.RTY2WDOT(q.rho, q.T, q.Y, wdot);
  amrex::Real wmax = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    if (!amrex::Math::isfinite(wdot[n])) {
      ok = false;
      return 0.0;
    }
    wmax = amrex::max<amrex::Real>(wmax, std::abs(wdot[n]));
  }
  if (wmax <= 0.0) {
    return 0.0; // frozen chemistry here; nothing to correct
  }

#if defined(USE_SRK_EOS) || defined(USE_MANIFOLD_EOS)
  // Directional finite difference along the reaction path.
  const amrex::Real dY_scale = wmax / q.rho; // [1/s]
  const amrex::Real tau = 1.0e-6 / dY_scale; // gives max |dY| ~ 1e-6
  amrex::Real Yp[NUM_SPECIES], sum = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    Yp[n] = amrex::max<amrex::Real>(q.Y[n] + tau * wdot[n] / q.rho, 0.0);
    sum += Yp[n];
  }
  if (!(sum > 0.0)) {
    ok = false;
    return 0.0;
  }
  for (int n = 0; n < NUM_SPECIES; n++) {
    Yp[n] /= sum;
  }
  amrex::Real Tp = q.T, pp = 0.0;
  eos.REY2T(q.rho, q.e, Yp, Tp);
  eos.RTY2P(q.rho, Tp, Yp, pp);
  if (!amrex::Math::isfinite(pp) || !amrex::Math::isfinite(Tp) || !(Tp > 0.0)) {
    ok = false;
    return 0.0;
  }
  return (pp - q.p) / tau;
#else
  amrex::Real imw[NUM_SPECIES], ei[NUM_SPECIES], cv = 0.0;
  eos.inv_molecular_weight(imw);
  eos.RTY2Ei(q.rho, q.T, q.Y, ei);
  eos.RTY2Cv(q.rho, q.T, q.Y, cv);
  if (!amrex::Math::isfinite(cv) || !(cv > 0.0)) {
    ok = false;
    return 0.0;
  }
  amrex::Real mole_term = 0.0, heat_term = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    mole_term += wdot[n] * imw[n];
    heat_term += ei[n] * wdot[n];
  }
  const amrex::Real dpdt = pele::physics::Constants::RU * q.T * mole_term -
                           q.p * heat_term / (q.rho * q.T * cv);
  if (!amrex::Math::isfinite(dpdt)) {
    ok = false;
    return 0.0;
  }
  return dpdt;
#endif
}

// ---------------------------------------------------------------------------
//  Pack a set of ghost primitives into a conserved state.
//
//  Every NVAR slot is written, and the algebraic identities that the rest of
//  the code assumes are enforced exactly rather than approximately:
//      sum_k (rho Y)_k == rho          (species renormalised before packing)
//      UEDEN == UEINT + 0.5 rho |u|^2
//  The passive families (advected, auxiliary, linear/soot) ride the lambda_0
//  characteristic and are handled by the caller through `pass`.
// ---------------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
pack_ghost(
  const amrex::Real rho,
  const amrex::Real u[3],
  const amrex::Real e,
  const amrex::Real T,
  const amrex::Real Y[NUM_SPECIES],
  const amrex::Real s_N[NVAR],
  amrex::Real s_ghost[NVAR]) noexcept
{
  s_ghost[URHO] = rho;
  s_ghost[UMX] = rho * u[0];
  s_ghost[UMY] = rho * u[1];
  s_ghost[UMZ] = rho * u[2];
  s_ghost[UEINT] = rho * e;
  s_ghost[UEDEN] =
    rho * e + 0.5 * rho * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
  s_ghost[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    s_ghost[UFS + n] = rho * Y[n];
  }
  // Passive families default to a mass-scaled zero-gradient; the caller
  // overwrites them where a target applies.
  const amrex::Real ratio =
    (s_N[URHO] > constants::very_small_num()) ? rho / s_N[URHO] : 1.0;
  amrex::ignore_unused(ratio);
#if NUM_ADV > 0
  for (int n = 0; n < NUM_ADV; n++) {
    s_ghost[UFA + n] = s_N[UFA + n] * ratio; // per unit mass in ctoprim
  }
#endif
#if NUM_AUX > 0
  for (int n = 0; n < NUM_AUX; n++) {
    s_ghost[UFX + n] = s_N[UFX + n]; // NOT per unit mass in ctoprim
  }
#endif
#if NUM_LIN > 0
  for (int n = 0; n < NUM_LIN; n++) {
    s_ghost[ULIN + n] = s_N[ULIN + n]; // NOT per unit mass in ctoprim
  }
#endif
}

// ===========================================================================
//  apply()
//
//  Fill one ghost cell.
//
//  s_N, s_Nm1, s_Nm2  interior stencil walking INWARD from the boundary.  The
//                     caller is responsible for clamping the TANGENTIAL
//                     components of these indices into the domain, so that the
//                     fill is a pure function of valid interior data and is
//                     therefore deterministic and free of any dependence on
//                     the order in which ghost cells are visited.
//  n_stencil          how many of the three are usable (3, 2 or 1).  Reduced
//                     by the caller, and further reduced here if a stencil
//                     cell turns out to hold EB body state.
//  dx_normal          cell size along the boundary normal.
//  idir               boundary-normal direction, 0/1/2.
//  sgn                +1 on a low-side face, -1 on a high-side face (AMReX
//                     convention, matching bcnormal).
//  layer              1 for the ghost cell nearest the domain, 2 for the next,
//                     and so on.
//  diag               optional amrex::Long counter array of length
//                     Diag::count; may be
//                     null.
// ===========================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
apply(
  const amrex::Real s_N[NVAR],
  const amrex::Real s_Nm1[NVAR],
  const amrex::Real s_Nm2[NVAR],
  int n_stencil,
  const amrex::Real dx_normal,
  const int idir,
  const int sgn,
  const int layer,
  Target const& tgt,
  Params const& prm,
  amrex::Real s_ghost[NVAR],
  amrex::Long* diag = nullptr,
  Transverse const* tr = nullptr) noexcept
{
  auto eos = pele::physics::PhysicsType::eos();

  auto bump = [=](const int which) {
    if (diag != nullptr) {
      amrex::HostDevice::Atomic::Add(&diag[which], amrex::Long(1));
    }
  };

  // ---- 0. Stencil hygiene ------------------------------------------------
  if (is_body_state(s_N)) {
    // Nothing usable to build a characteristic state from.  Leave whatever the
    // caller pre-populated (the reflected NoSlipWall state) in place.
    bump(Diag::body_state);
    copy_state(s_N, s_ghost);
    return;
  }
  if ((n_stencil > 2) && is_body_state(s_Nm2)) {
    n_stencil = 2;
    bump(Diag::body_state);
  }
  if ((n_stencil > 1) && is_body_state(s_Nm1)) {
    n_stencil = 1;
    bump(Diag::body_state);
  }
  // Extrapolation order can never exceed what the stencil supports.  A
  // minmod-limited slope needs three cells; with two it would be unlimited.
  const int order =
    amrex::min<int>(prm.order, amrex::max<int>(n_stencil - 1, 1));
  const bool linear = (order >= 2) && (n_stencil >= 3);

  const CellPrim qN = cell_primitives(s_N);
  if (!qN.ok) {
    bump(Diag::eos_failure);
    copy_state(s_N, s_ghost);
    return;
  }

  // ---- 1. Outward-normal frame ------------------------------------------
  const amrex::Real n_sgn = static_cast<amrex::Real>(-sgn);
  const amrex::Real u_out_N = n_sgn * qN.u[idir];
  const amrex::Real c_N = qN.c;
  const amrex::Real mach = u_out_N / c_N;

  // Frozen acoustic impedance (see the header note: this must be a single
  // value shared by every stencil cell).
  const amrex::Real rho_c = qN.rho * c_N;

  // ---- 2. Regimes with an exact answer ----------------------------------
  if (std::abs(mach) >= 1.0) {
    if (mach >= 1.0) {
      // Supersonic outflow: no characteristic can carry information inward.
      // Zero-gradient is exact, not an approximation.
      bump(Diag::supersonic);
      copy_state(s_N, s_ghost);
      return;
    }
    // Supersonic inflow: every characteristic is incoming, so the full state
    // is specified.  Only meaningful with an inflow target.
    if (tgt.type == Type::inflow) {
      amrex::Real rho_g = 0.0, e_g = 0.0;
      eos.PYT2RE(tgt.p > 0.0 ? tgt.p : qN.p, tgt.Y, tgt.T, rho_g, e_g);
      if (
        amrex::Math::isfinite(rho_g) && (rho_g > 0.0) &&
        amrex::Math::isfinite(e_g)) {
        amrex::Real u_g[3] = {tgt.u[0], tgt.u[1], tgt.u[2]};
        pack_ghost(rho_g, u_g, e_g, tgt.T, tgt.Y, s_N, s_ghost);
        bump(Diag::supersonic);
        return;
      }
      bump(Diag::eos_failure);
    }
    bump(Diag::supersonic);
    copy_state(s_N, s_ghost);
    return;
  }

  // ---- 3. Direction sanity ----------------------------------------------
  const bool outflow_face = (tgt.type == Type::outflow);
  if (outflow_face && (u_out_N < 0.0)) {
    // Transient reversal at an outflow.  Do not fight the flow with an
    // outflow model that assumes the lambda_0 family is leaving; fall back to
    // an ambient-pressure closure that keeps the interior composition and
    // tangential velocity so a grazing flow is not arrested.
    bump(Diag::reversed);
    amrex::Real T_g = qN.T, e_g = 0.0;
    eos.RYP2T(qN.rho, qN.Y, tgt.p, T_g);
    eos.RTY2E(qN.rho, T_g, qN.Y, e_g);
    if (
      !amrex::Math::isfinite(T_g) || !(T_g > 0.0) ||
      !amrex::Math::isfinite(e_g)) {
      bump(Diag::eos_failure);
      copy_state(s_N, s_ghost);
      return;
    }
    pack_ghost(qN.rho, qN.u, e_g, T_g, qN.Y, s_N, s_ghost);
    return;
  }
  if (!outflow_face && (u_out_N > 0.0)) {
    // Outflow through a face configured as an inflow.  The inflow model would
    // impose D+NUM_SPECIES+1 conditions on a point that admits one; that is
    // ill-posed.  Zero-gradient and count it.
    bump(Diag::reversed);
    copy_state(s_N, s_ghost);
    return;
  }

  // ---- 4. Outgoing acoustic invariant, R_+ ------------------------------
  // R_+ leaves the domain at both subsonic inflow and subsonic outflow
  // (lambda_+ = u_out + c > 0 in both cases), so it is always extrapolated.
  const amrex::Real Rp_N = u_out_N + qN.p / rho_c;
  amrex::Real dRp = 0.0;
  CellPrim qNm1, qNm2;
  if (linear) {
    qNm1 = cell_primitives(s_Nm1);
    qNm2 = cell_primitives(s_Nm2);
    if (qNm1.ok && qNm2.ok) {
      const amrex::Real Rp_Nm1 = n_sgn * qNm1.u[idir] + qNm1.p / rho_c;
      const amrex::Real Rp_Nm2 = n_sgn * qNm2.u[idir] + qNm2.p / rho_c;
      dRp = minmod(Rp_N - Rp_Nm1, Rp_Nm1 - Rp_Nm2);
    } else {
      bump(Diag::eos_failure);
    }
  }
  const amrex::Real Rp_g = Rp_N + static_cast<amrex::Real>(layer) * dRp;

  // ---- 5. Incoming acoustic invariant, R_- ------------------------------
  // lambda_- = u_out - c < 0 always (subsonic), so R_- is always modelled.
  //
  //     L_in       = K (phi - phi_target)
  //     dR_-/dn    = -L_in / (lambda_- rho c) = L_in / ((c - u_out) rho c)
  //     R_-,ghost  = R_-,N + layer dx dR_-/dn
  //
  // Sign check, outflow: phi = p, K > 0.  p_N > p_target gives dR_-/dn > 0,
  // hence R_-,ghost > R_-,N, and since dp/dR_- = -rho c / 2 the ghost pressure
  // falls: the boundary is pulled toward the target.  Anchoring, as intended.
  //
  // Sign check, inflow: phi = u_out with a leading minus (see below).  u_out_N
  // above the target gives dR_-/dn < 0, hence a smaller R_-,ghost, and since
  // du_out/dR_- = +1/2 the ghost normal velocity falls toward the target.
  const amrex::Real Rm_N = u_out_N - qN.p / rho_c;
  const amrex::Real lambda_in = c_N - u_out_N; // > 0, magnitude of the speed
  const amrex::Real one_m_M2 =
    amrex::max<amrex::Real>(1.0 - mach * mach, constants::very_small_num());

  // Transverse contribution to the modelled incoming wave.  Additive on the
  // amplitude, which is exactly why the relaxation is written as a gradient
  // increment: there is a well-defined slot for it here, and none at all in a
  // hard far-field pin of R_-.
  amrex::Real beta = prm.beta;
  if (beta < 0.0) {
    beta = amrex::min<amrex::Real>(std::abs(mach), 1.0);
  }
  amrex::Real T_in = 0.0;
  if ((tr != nullptr) && tr->valid && (beta < 1.0)) {
    amrex::Real gam = 0.0;
    eos.RTY2G(qN.rho, qN.T, qN.Y, gam);
    amrex::Real div_ut = 0.0;
    bool any = false;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      if (
        (d == idir) || (tr->sm[d] == nullptr) || (tr->sp[d] == nullptr) ||
        (tr->inv_dt[d] == 0.0)) {
        continue;
      }
      const CellPrim qm = cell_primitives(tr->sm[d]);
      const CellPrim qp = cell_primitives(tr->sp[d]);
      if (!qm.ok || !qp.ok) {
        bump(Diag::transverse_drop);
        continue;
      }
      const amrex::Real dpdt = (qp.p - qm.p) * tr->inv_dt[d];
      // du_out/dt: the normal velocity in the OUTWARD frame, differentiated
      // along the tangential direction.  The n_sgn here is what turns the
      // legacy hi-face T1 and lo-face T4 into one expression.
      const amrex::Real duodt =
        n_sgn * (qp.u[idir] - qm.u[idir]) * tr->inv_dt[d];
      const amrex::Real dutdt = (qp.u[d] - qm.u[d]) * tr->inv_dt[d];
      T_in += qN.u[d] * (dpdt - rho_c * duodt);
      div_ut += dutdt;
      any = true;
    }
    if (any && amrex::Math::isfinite(gam) && (gam > 0.0)) {
      T_in += gam * qN.p * div_ut;
    } else {
      T_in = 0.0;
      bump(Diag::transverse_drop);
    }
    if (!amrex::Math::isfinite(T_in)) {
      T_in = 0.0;
      bump(Diag::transverse_drop);
    }
  }

  // Entropy-family slope at the face: measured for two purposes below.  dS
  // carries no acoustic content and is an OUTGOING family at an outflow, so
  // nothing the incoming model consumes here can feed on the waves it
  // launches.  Measured on the same stencil, in the same frozen-impedance
  // variables, and with the same limiter as R_+.
  amrex::Real dRm = 0.0;
  if (outflow_face && linear && qNm1.ok && qNm2.ok) {
    const amrex::Real inv_c2 = 1.0 / (c_N * c_N);
    const amrex::Real S_N = qN.rho - qN.p * inv_c2;
    const amrex::Real S_Nm1 = qNm1.rho - qNm1.p * inv_c2;
    const amrex::Real S_Nm2 = qNm2.rho - qNm2.p * inv_c2;
    const amrex::Real dS = minmod(S_N - S_Nm1, S_Nm1 - S_Nm2);

    // The transit guard.  Material structure sitting in the outflow boundary
    // cell is the configuration for which no dial at its default survives a
    // crossing (the flame-exit table in NSCBC-FlameOutflow/README.md): count
    // it and let nscbc_report_diagnostics() say so.  5% of rho per cell is
    // well above acoustic or roundoff content and well below any front.
    if ((layer == 1) && (std::abs(dS) > 0.05 * qN.rho)) {
      bump(Diag::structure);
    }

    // Material-slope continuation of R_- (see Params::extrap_material): the
    // measured slope of R_- is used only within the bound the entropy family
    // supplies through steady continuity and the momentum-flux pressure
    // gradient.
    if (prm.extrap_material) {
      const amrex::Real Rm_Nm1 = n_sgn * qNm1.u[idir] - qNm1.p / rho_c;
      const amrex::Real Rm_Nm2 = n_sgn * qNm2.u[idir] - qNm2.p / rho_c;
      const amrex::Real dRm_meas = minmod(Rm_N - Rm_Nm1, Rm_Nm1 - Rm_Nm2);
      const amrex::Real du_mat = -u_out_N * dS / (qN.rho * one_m_M2);
      dRm = minmod(dRm_meas, du_mat * (1.0 + mach));
      if (!amrex::Math::isfinite(dRm)) {
        dRm = 0.0;
        bump(Diag::floored);
      }
    }
  }

  amrex::Real Rm_g;
  if (outflow_face && prm.pin_farfield) {
    // Hard far-field pin: R_- takes its far-field value (u_out = 0,
    // p = p_target).  Value constraint, not a rate, and there is nowhere in it
    // to put a transverse term -- so pin_farfield ignores beta.
    Rm_g = -tgt.p / rho_c;
  } else {
    amrex::Real L_in;
    if (outflow_face) {
      const amrex::Real K = prm.sigma * one_m_M2 * c_N / prm.L_ref;
      L_in = K * (qN.p - tgt.p);
    } else {
      const amrex::Real u_out_tgt = n_sgn * tgt.u[idir];
      const amrex::Real K = prm.relax_u * one_m_M2 * c_N / prm.L_ref;
      L_in = -K * qN.rho * c_N * (u_out_N - u_out_tgt);
    }
    L_in -= (1.0 - beta) * T_in;
    if (prm.beta_s < 1.0) {
      bool src_ok = true;
      const amrex::Real S_p = reaction_dpdt(qN, src_ok);
      if (src_ok) {
        L_in -= (1.0 - prm.beta_s) * S_p;
      } else {
        bump(Diag::source_drop);
      }
    }
    // There is deliberately NO diffusive source term here, and that is a
    // measured decision, not an omission.  A dp/dt|_diffusion of the boundary
    // cell -- normal conduction and mixture-averaged species diffusion of the
    // resolved fields, converted through the EOS -- was built, verified exact
    // on quadratic profiles, and measured: with real conduction in the 1-D
    // driver the boundary error is 887 dyn/cm2 under the entropy closure,
    // 104 under extrap_temperature, and adding the modelled term moves 104 to
    // -911; in PeleC's flame-outflow at sigma = 1 it moves +1200 to +1771.
    // The ghost-cell form does not HAVE the flux-form's viscous-condition
    // gap: the diffusion operator reads these ghost cells, so a correct
    // ghost closure (extrap_temperature, gated by C8 and C12) already
    // carries the diffusive physics, and an amplitude-side term counts it
    // twice.  See check C12 in Verification/NSCBC1D.
    Rm_g = Rm_N + static_cast<amrex::Real>(layer) * dRm +
           static_cast<amrex::Real>(layer) * dx_normal * L_in /
             (lambda_in * rho_c);
  }

  // ---- 6. Invert the acoustic pair --------------------------------------
  amrex::Real u_out_g = 0.5 * (Rp_g + Rm_g);
  amrex::Real p_g = 0.5 * rho_c * (Rp_g - Rm_g);

  const amrex::Real p_floor = amrex::max<amrex::Real>(
    1.0e-3 * (tgt.p > 0.0 ? tgt.p : qN.p), constants::very_small_num());
  if (!amrex::Math::isfinite(p_g) || (p_g < p_floor)) {
    p_g = p_floor;
    bump(Diag::floored);
  }
  if (!amrex::Math::isfinite(u_out_g)) {
    u_out_g = u_out_N;
    bump(Diag::floored);
  }

  amrex::Real u_g[3] = {qN.u[0], qN.u[1], qN.u[2]};
  u_g[idir] = n_sgn * u_out_g;

  // ---- 7. The lambda_0 family -------------------------------------------
  // Entropy, tangential velocities, species and passive scalars all convect at
  // u_out.  At an OUTFLOW they leave the domain and are extrapolated; at an
  // INFLOW they enter and are set from the target -- composition hard (all
  // NUM_SPECIES characteristics are incoming, so it is well posed), and
  // temperature and tangential velocity through a relaxation.
  amrex::Real rho_g = 0.0, T_g = 0.0, e_g = 0.0;
  amrex::Real Y_g[NUM_SPECIES];

  if (outflow_face) {
    // Entropy invariant S = rho - p/c^2, extrapolated with the same limited
    // linear form as R_+.  With a zero slope this reduces exactly to the
    // linearised isentrope rho_g = rho_N + (p_g - p_N)/c^2.
    const amrex::Real inv_c2 = 1.0 / (c_N * c_N);
    const amrex::Real S_N = qN.rho - qN.p * inv_c2;
    amrex::Real dS = 0.0;
    amrex::Real dY[NUM_SPECIES] = {0.0};
    amrex::Real dut[3] = {0.0, 0.0, 0.0};
    if (linear && qNm1.ok && qNm2.ok) {
      const amrex::Real S_Nm1 = qNm1.rho - qNm1.p * inv_c2;
      const amrex::Real S_Nm2 = qNm2.rho - qNm2.p * inv_c2;
      dS = minmod(S_N - S_Nm1, S_Nm1 - S_Nm2);
      for (int n = 0; n < NUM_SPECIES; n++) {
        dY[n] = minmod(qN.Y[n] - qNm1.Y[n], qNm1.Y[n] - qNm2.Y[n]);
      }
      for (int d = 0; d < 3; d++) {
        if (d != idir) {
          dut[d] = minmod(qN.u[d] - qNm1.u[d], qNm1.u[d] - qNm2.u[d]);
        }
      }
    }
    const amrex::Real fl = static_cast<amrex::Real>(layer);
    const amrex::Real S_g = S_N + fl * dS;
    rho_g = S_g + p_g * inv_c2;

    // Temperature closure: extrapolate T on the same limited slope and let rho
    // follow from the EOS below, so that the face temperature gradient is the
    // interior one rather than a by-product of the acoustic algebra.
    amrex::Real dT = 0.0;
    if (prm.extrap_temperature && linear && qNm1.ok && qNm2.ok) {
      dT = minmod(qN.T - qNm1.T, qNm1.T - qNm2.T);
    }

    amrex::Real Ysum = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      Y_g[n] = amrex::max<amrex::Real>(qN.Y[n] + fl * dY[n], 0.0);
      Ysum += Y_g[n];
    }
    if (Ysum > constants::very_small_num()) {
      const amrex::Real inv = 1.0 / Ysum;
      for (int n = 0; n < NUM_SPECIES; n++) {
        Y_g[n] *= inv;
      }
    } else {
      for (int n = 0; n < NUM_SPECIES; n++) {
        Y_g[n] = qN.Y[n];
      }
      bump(Diag::floored);
    }

    for (int d = 0; d < 3; d++) {
      if (d != idir) {
        u_g[d] = qN.u[d] + fl * dut[d];
      }
    }

    if (prm.extrap_temperature) {
      T_g = qN.T + fl * dT;
      // A temperature floor rather than a density one: rho is the derived
      // quantity in this closure.
      const amrex::Real T_floor =
        amrex::max<amrex::Real>(1.0e-3 * qN.T, constants::very_small_num());
      if (!amrex::Math::isfinite(T_g) || (T_g < T_floor)) {
        T_g = T_floor;
        bump(Diag::floored);
      }
      eos.PYT2RE(p_g, Y_g, T_g, rho_g, e_g);
    } else {
      const amrex::Real rho_floor =
        amrex::max<amrex::Real>(1.0e-6 * qN.rho, constants::very_small_num());
      if (!amrex::Math::isfinite(rho_g) || (rho_g < rho_floor)) {
        rho_g = rho_floor;
        bump(Diag::floored);
      }
      eos.RYP2T(rho_g, Y_g, p_g, T_g);
      eos.RTY2E(rho_g, T_g, Y_g, e_g);
    }
  } else {
    // Inflow.  A single dimensionless nudge factor is shared by the
    // temperature and the tangential velocities: both ride lambda_0, whose
    // magnitude is |u_out|.
    //
    //     nu = layer dx relax_t c / (L_ref |u_out|)   clamped to [0, 1]
    //
    // The legacy Fortran used ratios rho c R_specific / L for temperature and
    // c / L for tangential velocity, which are not dimensionally comparable
    // and are why relax_T had to be given a negative value.  One consistent
    // factor, applied as an explicit move toward the target, cannot get the
    // sign wrong.
    const amrex::Real u_mag =
      amrex::max<amrex::Real>(std::abs(u_out_N), constants::smallu());
    const amrex::Real nu = amrex::min<amrex::Real>(
      static_cast<amrex::Real>(layer) * dx_normal * prm.relax_t * c_N /
        (prm.L_ref * u_mag),
      1.0);

    T_g = qN.T - nu * (qN.T - tgt.T);
    for (int d = 0; d < 3; d++) {
      if (d != idir) {
        u_g[d] = qN.u[d] - nu * (qN.u[d] - tgt.u[d]);
      }
    }
    // Composition: hard.  Every species characteristic is incoming.
    amrex::Real Ysum = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      Y_g[n] = amrex::max<amrex::Real>(tgt.Y[n], 0.0);
      Ysum += Y_g[n];
    }
    if (Ysum > constants::very_small_num()) {
      const amrex::Real inv = 1.0 / Ysum;
      for (int n = 0; n < NUM_SPECIES; n++) {
        Y_g[n] *= inv;
      }
    } else {
      for (int n = 0; n < NUM_SPECIES; n++) {
        Y_g[n] = qN.Y[n];
      }
      bump(Diag::floored);
    }
    if (!amrex::Math::isfinite(T_g) || !(T_g > 0.0)) {
      T_g = qN.T;
      bump(Diag::floored);
    }
    // Close on (p, T, Y): p from the acoustic pair, T and Y from the target.
    eos.PYT2RE(p_g, Y_g, T_g, rho_g, e_g);
  }

  if (
    !amrex::Math::isfinite(rho_g) || !(rho_g > 0.0) ||
    !amrex::Math::isfinite(e_g) || !amrex::Math::isfinite(T_g) ||
    !(T_g > 0.0)) {
    bump(Diag::eos_failure);
    copy_state(s_N, s_ghost);
    return;
  }

  // ---- 8. Pack and sanitise ---------------------------------------------
  pack_ghost(rho_g, u_g, e_g, T_g, Y_g, s_N, s_ghost);

  for (int n = 0; n < NVAR; n++) {
    if (!amrex::Math::isfinite(s_ghost[n])) {
      s_ghost[n] = s_N[n];
      bump(Diag::floored);
    }
  }
}

} // namespace pc_nscbc

#endif
