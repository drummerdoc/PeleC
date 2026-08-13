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
//  Outgoing quantities (R_+ always; S, u_t, Y, scalars at an outflow) are
//  extrapolated outward from the interior with a minmod-limited linear slope.
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
//    1. Transverse terms ARE available, weighted by Params::beta -- see the
//       guidance on that member; the optimum is near 0.5 and the failure mode
//       below ~0.3 is severe.  The reaction and viscous source contributions
//       are NOT available: they need sources_for_hydro, which does not exist
//       during FillPatch.  With
//       beta = 1 (the default) the transverse terms are discarded too and this
//       is the pure 1-D-normal LODI limit.  Either way, place an outflow away
//       from flames: heat release in the boundary cell has no model here and
//       shows up as a mean pressure offset.
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
  // of R_- to the far-field value (u_out = 0, p = p_target).  This is the
  // formulation used by CAMR's PS_nscbc.H.  It is simultaneously non-reflecting
  // and anchored, but it anchors to p_target + rho c u_out rather than to
  // p_target, and it is a value constraint rather than a rate, so its effective
  // relaxation rate is c/dx -- i.e. it does not converge under grid refinement.
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
  amrex::Real p = 0.0;                // outflow: target pressure [Pa]
  amrex::Real u[3] = {0.0, 0.0, 0.0}; // inflow: target velocity [m/s]
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
//  array reduced at the end of the run.
// ---------------------------------------------------------------------------
namespace Diag {
enum : int {
  supersonic = 0,  // |M| >= 1: zero-gradient copy is exact, not a failure
  reversed,        // flow ran the wrong way through the face
  body_state,      // an EB-covered cell was found in the stencil
  eos_failure,     // an EOS call returned a non-finite or non-positive value
  floored,         // a positivity floor was hit
  transverse_drop, // a transverse derivative was unavailable or unusable
  count
};
}

// ===========================================================================
//  Small helpers
// ===========================================================================

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
finite_or(const amrex::Real x, const amrex::Real fallback) noexcept
{
  return std::isfinite(x) ? x : fallback;
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
  return !(s[URHO] > constants::very_small_num()) || !std::isfinite(s[URHO]) ||
         !std::isfinite(s[UEDEN]);
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
  if (!(q.T > 0.0) || !std::isfinite(q.T)) {
    q.T = 300.0; // initial guess for the iterative fallback
    eos.REY2T(q.rho, q.e, q.Y, q.T);
  }

  eos.RTY2P(q.rho, q.T, q.Y, q.p);
  eos.RTY2Cs(q.rho, q.T, q.Y, q.c);

  if (
    !std::isfinite(q.p) || !(q.p > 0.0) || !std::isfinite(q.c) ||
    !(q.c > 0.0) || !std::isfinite(q.T) || !(q.T > 0.0)) {
    q.ok = false;
  }
  return q;
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
//  diag               optional counter array of length Diag::count; may be
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
  int* diag = nullptr,
  Transverse const* tr = nullptr) noexcept
{
  auto eos = pele::physics::PhysicsType::eos();

  auto bump = [=](const int which) {
    if (diag != nullptr) {
      amrex::HostDevice::Atomic::Add(&diag[which], 1);
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
      if (std::isfinite(rho_g) && (rho_g > 0.0) && std::isfinite(e_g)) {
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
    if (!std::isfinite(T_g) || !(T_g > 0.0) || !std::isfinite(e_g)) {
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
    if (any && std::isfinite(gam) && (gam > 0.0)) {
      T_in += gam * qN.p * div_ut;
    } else {
      T_in = 0.0;
      bump(Diag::transverse_drop);
    }
    if (!std::isfinite(T_in)) {
      T_in = 0.0;
      bump(Diag::transverse_drop);
    }
  }

  amrex::Real Rm_g;
  if (outflow_face && prm.pin_farfield) {
    // CAMR-style hard pin: R_- takes its far-field value (u_out = 0,
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
    Rm_g = Rm_N + static_cast<amrex::Real>(layer) * dx_normal * L_in /
                    (lambda_in * rho_c);
  }

  // ---- 6. Invert the acoustic pair --------------------------------------
  amrex::Real u_out_g = 0.5 * (Rp_g + Rm_g);
  amrex::Real p_g = 0.5 * rho_c * (Rp_g - Rm_g);

  const amrex::Real p_floor = amrex::max<amrex::Real>(
    1.0e-3 * (tgt.p > 0.0 ? tgt.p : qN.p), constants::very_small_num());
  if (!std::isfinite(p_g) || (p_g < p_floor)) {
    p_g = p_floor;
    bump(Diag::floored);
  }
  if (!std::isfinite(u_out_g)) {
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

    const amrex::Real rho_floor =
      amrex::max<amrex::Real>(1.0e-6 * qN.rho, constants::very_small_num());
    if (!std::isfinite(rho_g) || (rho_g < rho_floor)) {
      rho_g = rho_floor;
      bump(Diag::floored);
    }
    eos.RYP2T(rho_g, Y_g, p_g, T_g);
    eos.RTY2E(rho_g, T_g, Y_g, e_g);
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
    if (!std::isfinite(T_g) || !(T_g > 0.0)) {
      T_g = qN.T;
      bump(Diag::floored);
    }
    // Close on (p, T, Y): p from the acoustic pair, T and Y from the target.
    eos.PYT2RE(p_g, Y_g, T_g, rho_g, e_g);
  }

  if (
    !std::isfinite(rho_g) || !(rho_g > 0.0) || !std::isfinite(e_g) ||
    !std::isfinite(T_g) || !(T_g > 0.0)) {
    bump(Diag::eos_failure);
    copy_state(s_N, s_ghost);
    return;
  }

  // ---- 8. Pack and sanitise ---------------------------------------------
  pack_ghost(rho_g, u_g, e_g, T_g, Y_g, s_N, s_ghost);

  for (int n = 0; n < NVAR; n++) {
    if (!std::isfinite(s_ghost[n])) {
      s_ghost[n] = s_N[n];
      bump(Diag::floored);
    }
  }
}

} // namespace pc_nscbc

#endif
