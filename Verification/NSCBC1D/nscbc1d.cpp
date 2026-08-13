// ============================================================================
//  nscbc1d -- standalone 1-D verification driver for Source/NSCBC.H
//
//  This driver compiles the production NSCBC kernel unmodified against a
//  minimal 1-D finite-volume Euler solver.  Its purpose is to separate two
//  failure modes that are otherwise indistinguishable in a full PeleC run:
//
//      "the boundary condition is wrong"   vs   "the AMReX plumbing is wrong"
//
//  It is cheap enough to sweep parameters, so it also produces the reference
//  curves (reflection coefficient vs sigma, relaxation rate vs resolution)
//  that the AMReX-side regression tests are checked against.
//
//  Build:  see CMakeLists.txt in this directory.
//  Run:    ./nscbc1d              (runs every check, prints a summary)
//          ./nscbc1d sweep        (also dumps the sigma sweep as CSV)
//
//  Solver: MUSCL-Hancock, minmod-limited primitive reconstruction, HLLC flux,
//  SSP-RK2.  Deliberately simple -- it is the boundary condition under test,
//  not the interior scheme.
// ============================================================================

#include <AMReX.H>
#include <AMReX_Print.H>

#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#include "NSCBC.H"

using amrex::Real;

namespace {

constexpr int NG = 4; // ghost layers, matching PeleC::numGrow() without EB

// ---------------------------------------------------------------------------
//  A 1-D state vector, laid out exactly like PeleC's conserved state so that
//  the kernel sees the indices it expects.
// ---------------------------------------------------------------------------
struct Field
{
  int n;
  std::vector<Real> s; // (n + 2*NG) * NVAR

  explicit Field(int n_) : n(n_), s(static_cast<size_t>(n + 2 * NG) * NVAR, 0.0)
  {
  }
  Real* at(int i) { return &s[static_cast<size_t>(i + NG) * NVAR]; }
  const Real* at(int i) const { return &s[static_cast<size_t>(i + NG) * NVAR]; }
};

struct Prim
{
  Real rho, u, p, T;
  Real Y[NUM_SPECIES];
};

void
set_state(Real* s, Real rho, Real u, Real T, const Real Y[NUM_SPECIES])
{
  auto eos = pele::physics::PhysicsType::eos();
  Real e = 0.0;
  eos.RTY2E(rho, T, Y, e);
  s[URHO] = rho;
  s[UMX] = rho * u;
  s[UMY] = 0.0;
  s[UMZ] = 0.0;
  s[UEINT] = rho * e;
  s[UEDEN] = rho * (e + 0.5 * u * u);
  s[UTEMP] = T;
  for (int k = 0; k < NUM_SPECIES; k++) {
    s[UFS + k] = rho * Y[k];
  }
}

Prim
get_prim(const Real* s)
{
  auto eos = pele::physics::PhysicsType::eos();
  Prim q{};
  q.rho = s[URHO];
  q.u = s[UMX] / q.rho;
  q.T = s[UTEMP];
  Real ysum = 0.0;
  for (int k = 0; k < NUM_SPECIES; k++) {
    q.Y[k] = s[UFS + k] / q.rho;
    ysum += q.Y[k];
  }
  for (int k = 0; k < NUM_SPECIES; k++) {
    q.Y[k] /= ysum;
  }
  eos.RTY2P(q.rho, q.T, q.Y, q.p);
  return q;
}

Real
sound_speed(const Prim& q)
{
  auto eos = pele::physics::PhysicsType::eos();
  Real c = 0.0;
  eos.RTY2Cs(q.rho, q.T, q.Y, c);
  return c;
}

// ---------------------------------------------------------------------------
//  Boundary fill: exactly the dispatch the AMReX-side BCfill.cpp will do.
//  layer runs 1..NG outward; the stencil walks inward from the boundary cell.
// ---------------------------------------------------------------------------
void
fill_bcs(
  Field& f,
  const pc_nscbc::Target& lo_tgt,
  const pc_nscbc::Target& hi_tgt,
  const pc_nscbc::Params& prm,
  Real dx,
  int* diag = nullptr)
{
  const int n = f.n;
  for (int layer = 1; layer <= NG; layer++) {
    if (lo_tgt.type != pc_nscbc::Type::off) {
      pc_nscbc::apply(
        f.at(0), f.at(1), f.at(2), 3, dx, /*idir=*/0, /*sgn=*/+1, layer, lo_tgt,
        prm, f.at(-layer), diag);
    } else {
      for (int v = 0; v < NVAR; v++) {
        f.at(-layer)[v] = f.at(0)[v];
      }
    }
    if (hi_tgt.type != pc_nscbc::Type::off) {
      pc_nscbc::apply(
        f.at(n - 1), f.at(n - 2), f.at(n - 3), 3, dx, /*idir=*/0, /*sgn=*/-1,
        layer, hi_tgt, prm, f.at(n - 1 + layer), diag);
    } else {
      for (int v = 0; v < NVAR; v++) {
        f.at(n - 1 + layer)[v] = f.at(n - 1)[v];
      }
    }
  }
}

// ---------------------------------------------------------------------------
//  HLLC flux for the multi-species Euler equations.
// ---------------------------------------------------------------------------
void
hllc(const Prim& L, const Prim& R, Real* flx)
{
  auto eos = pele::physics::PhysicsType::eos();
  const Real cL = sound_speed(L), cR = sound_speed(R);
  const Real sL = std::min(L.u - cL, R.u - cR);
  const Real sR = std::max(L.u + cL, R.u + cR);
  const Real sM =
    (R.p - L.p + L.rho * L.u * (sL - L.u) - R.rho * R.u * (sR - R.u)) /
    (L.rho * (sL - L.u) - R.rho * (sR - R.u));

  auto cons_flux = [&](const Prim& q, Real* F) {
    Real e = 0.0;
    eos.RTY2E(q.rho, q.T, q.Y, e);
    const Real E = q.rho * (e + 0.5 * q.u * q.u);
    F[URHO] = q.rho * q.u;
    F[UMX] = q.rho * q.u * q.u + q.p;
    F[UMY] = 0.0;
    F[UMZ] = 0.0;
    F[UEDEN] = (E + q.p) * q.u;
    F[UEINT] = q.rho * e * q.u;
    F[UTEMP] = 0.0;
    for (int k = 0; k < NUM_SPECIES; k++) {
      F[UFS + k] = q.rho * q.Y[k] * q.u;
    }
  };
  auto star_state = [&](const Prim& q, Real s, Real* U) {
    Real e = 0.0;
    eos.RTY2E(q.rho, q.T, q.Y, e);
    const Real E = q.rho * (e + 0.5 * q.u * q.u);
    const Real fac = (s - q.u) / (s - sM);
    const Real rs = q.rho * fac;
    U[URHO] = rs;
    U[UMX] = rs * sM;
    U[UMY] = 0.0;
    U[UMZ] = 0.0;
    U[UEDEN] = rs * (E / q.rho + (sM - q.u) * (sM + q.p / (q.rho * (s - q.u))));
    U[UEINT] = rs * e;
    U[UTEMP] = 0.0;
    for (int k = 0; k < NUM_SPECIES; k++) {
      U[UFS + k] = rs * q.Y[k];
    }
  };

  Real FL[NVAR], FR[NVAR], UL[NVAR], UR[NVAR], US[NVAR];
  cons_flux(L, FL);
  cons_flux(R, FR);
  auto cons_state = [&](const Prim& q, Real* U) {
    Real e = 0.0;
    eos.RTY2E(q.rho, q.T, q.Y, e);
    U[URHO] = q.rho;
    U[UMX] = q.rho * q.u;
    U[UMY] = 0.0;
    U[UMZ] = 0.0;
    U[UEDEN] = q.rho * (e + 0.5 * q.u * q.u);
    U[UEINT] = q.rho * e;
    U[UTEMP] = 0.0;
    for (int k = 0; k < NUM_SPECIES; k++) {
      U[UFS + k] = q.rho * q.Y[k];
    }
  };
  cons_state(L, UL);
  cons_state(R, UR);

  if (sL >= 0.0) {
    for (int v = 0; v < NVAR; v++) {
      flx[v] = FL[v];
    }
  } else if (sR <= 0.0) {
    for (int v = 0; v < NVAR; v++) {
      flx[v] = FR[v];
    }
  } else if (sM >= 0.0) {
    star_state(L, sL, US);
    for (int v = 0; v < NVAR; v++) {
      flx[v] = FL[v] + sL * (US[v] - UL[v]);
    }
  } else {
    star_state(R, sR, US);
    for (int v = 0; v < NVAR; v++) {
      flx[v] = FR[v] + sR * (US[v] - UR[v]);
    }
  }
}

Real
mm(Real a, Real b)
{
  return (a * b <= 0.0) ? 0.0 : ((std::abs(a) < std::abs(b)) ? a : b);
}

// One SSP-RK2 stage: conserved-variable MUSCL with minmod slopes.
void
stage(const Field& in, Field& out, Real dx, Real dt)
{
  const int n = in.n;
  std::vector<Real> flux(static_cast<size_t>(n + 1) * NVAR, 0.0);
  for (int i = 0; i <= n; i++) {
    Real sl[NVAR], sr[NVAR];
    for (int v = 0; v < NVAR; v++) {
      const Real dL =
        mm(in.at(i - 1)[v] - in.at(i - 2)[v], in.at(i)[v] - in.at(i - 1)[v]);
      const Real dR =
        mm(in.at(i)[v] - in.at(i - 1)[v], in.at(i + 1)[v] - in.at(i)[v]);
      sl[v] = in.at(i - 1)[v] + 0.5 * dL;
      sr[v] = in.at(i)[v] - 0.5 * dR;
    }
    // Recover T for each face state so the EOS is self-consistent.
    auto eos = pele::physics::PhysicsType::eos();
    auto to_prim = [&](Real* s) {
      Prim q{};
      q.rho = std::max(s[URHO], 1e-12);
      q.u = s[UMX] / q.rho;
      Real ys = 0.0;
      for (int k = 0; k < NUM_SPECIES; k++) {
        q.Y[k] = std::max(s[UFS + k] / q.rho, 0.0);
        ys += q.Y[k];
      }
      for (int k = 0; k < NUM_SPECIES; k++) {
        q.Y[k] /= ys;
      }
      const Real e = s[UEDEN] / q.rho - 0.5 * q.u * q.u;
      q.T = s[UTEMP] > 0.0 ? s[UTEMP] : 300.0;
      eos.REY2T(q.rho, e, q.Y, q.T);
      eos.RTY2P(q.rho, q.T, q.Y, q.p);
      return q;
    };
    hllc(to_prim(sl), to_prim(sr), &flux[static_cast<size_t>(i) * NVAR]);
  }
  for (int i = 0; i < n; i++) {
    for (int v = 0; v < NVAR; v++) {
      out.at(i)[v] =
        in.at(i)[v] - dt / dx *
                        (flux[static_cast<size_t>(i + 1) * NVAR + v] -
                         flux[static_cast<size_t>(i) * NVAR + v]);
    }
    // keep UTEMP consistent
    auto eos = pele::physics::PhysicsType::eos();
    Real* s = out.at(i);
    const Real rho = s[URHO];
    const Real u = s[UMX] / rho;
    Real Y[NUM_SPECIES], ys = 0.0;
    for (int k = 0; k < NUM_SPECIES; k++) {
      Y[k] = std::max(s[UFS + k] / rho, 0.0);
      ys += Y[k];
    }
    for (int k = 0; k < NUM_SPECIES; k++) {
      Y[k] /= ys;
    }
    const Real e = s[UEDEN] / rho - 0.5 * u * u;
    Real T = s[UTEMP] > 0.0 ? s[UTEMP] : 300.0;
    eos.REY2T(rho, e, Y, T);
    s[UTEMP] = T;
    s[UEINT] = rho * e;
  }
}

// PeleC and PelePhysics work in CGS: lengths in cm, pressures in dyn/cm^2,
// velocities in cm/s, densities in g/cm^3.  Constants::PATM = 1.01325e6.
struct Case
{
  int n = 400;
  Real L = 10.0;       // cm
  Real p0 = 1.01325e6; // dyn/cm^2 (1 atm)
  Real T0 = 300.0;     // K
  Real u0 = 0.0;       // cm/s
  Real cfl = 0.4;
};

// ---------------------------------------------------------------------------
//  Reporting
// ---------------------------------------------------------------------------
int n_pass = 0, n_fail = 0;

void
check(bool ok, const std::string& name, const std::string& detail)
{
  if (ok) {
    n_pass++;
    std::printf("  PASS  %-46s %s\n", name.c_str(), detail.c_str());
  } else {
    n_fail++;
    std::printf("  FAIL  %-46s %s\n", name.c_str(), detail.c_str());
  }
}

Real
air_Y(int k)
{
  return (k == O2_ID) ? 0.233 : 0.767;
}

} // namespace

// ===========================================================================
//  Checks
// ===========================================================================

// C1: a uniform state must be reproduced exactly in every ghost layer, for
//     both inflow and outflow, at any sigma.  If this fails, the kernel is
//     manufacturing a gradient out of nothing and everything downstream is
//     meaningless.
void
check_uniform()
{
  Case cs;
  Field f(cs.n);
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
  auto eos = pele::physics::PhysicsType::eos();
  Real rho0 = 0.0, e0 = 0.0;
  eos.PYT2RE(cs.p0, Y, cs.T0, rho0, e0);
  for (int i = -NG; i < cs.n + NG; i++) {
    set_state(f.at(i), rho0, cs.u0 + 2.0e3, cs.T0, Y); // 20 m/s
  }
  const Real dx = cs.L / cs.n;

  pc_nscbc::Params prm;
  prm.L_ref = cs.L;
  pc_nscbc::Target off, out, in;
  out.type = pc_nscbc::Type::outflow;
  out.p = cs.p0;
  in.type = pc_nscbc::Type::inflow;
  in.u[0] = 2.0e3;
  in.T = cs.T0;
  for (int k = 0; k < NUM_SPECIES; k++) {
    in.Y[k] = Y[k];
  }

  for (Real sig : {0.0, 0.25, 1.0}) {
    prm.sigma = sig;
    prm.relax_u = sig;
    prm.relax_t = sig;
    Field g = f;
    fill_bcs(g, in, out, prm, dx);
    Real worst = 0.0;
    for (int layer = 1; layer <= NG; layer++) {
      for (int v = 0; v < NVAR; v++) {
        const Real ref = std::max(std::abs(f.at(0)[v]), 1.0);
        worst = std::max(worst, std::abs(g.at(-layer)[v] - f.at(0)[v]) / ref);
        worst = std::max(
          worst, std::abs(g.at(cs.n - 1 + layer)[v] - f.at(cs.n - 1)[v]) / ref);
      }
    }
    char buf[160];
    std::snprintf(
      buf, sizeof(buf), "sigma=%.2f  max rel ghost error = %.3e", sig, worst);
    check(worst < 1e-12, "uniform state is reproduced exactly", buf);
  }
}

// C2: the relaxation must move the boundary TOWARD the target, for every
//     target quantity and on both faces.  This is the assertion that catches
//     the sign errors the legacy Fortran leaked to its users as "relax_T must
//     be negative".
void
check_relaxation_signs()
{
  Case cs;
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
  auto eos = pele::physics::PhysicsType::eos();
  const Real dx = cs.L / cs.n;
  pc_nscbc::Params prm;
  prm.L_ref = cs.L;
  prm.sigma = 0.25;
  prm.relax_u = 2.0;
  prm.relax_t = 0.2;
  prm.order = 1; // isolate the relaxation from the extrapolation

  // --- outflow: interior pressure above target -> ghost pressure below the
  //     non-relaxed value (i.e. pulled toward the target)
  for (int sgn : {+1, -1}) {
    Real rho = 0.0, e = 0.0;
    const Real p_int = 1.05 * cs.p0;
    eos.PYT2RE(p_int, Y, cs.T0, rho, e);
    Real sN[NVAR], sg_relaxed[NVAR], sg_free[NVAR];
    set_state(sN, rho, (sgn == +1) ? -3.0e3 : 3.0e3, cs.T0, Y); // outflow
    pc_nscbc::Target t;
    t.type = pc_nscbc::Type::outflow;
    t.p = cs.p0;
    pc_nscbc::apply(sN, sN, sN, 1, dx, 0, sgn, 1, t, prm, sg_relaxed);
    pc_nscbc::Params p0prm = prm;
    p0prm.sigma = 0.0;
    pc_nscbc::apply(sN, sN, sN, 1, dx, 0, sgn, 1, t, p0prm, sg_free);
    const Prim qr = get_prim(sg_relaxed), qf = get_prim(sg_free);
    char buf[160];
    std::snprintf(
      buf, sizeof(buf), "%s face: p_ghost %.3f vs %.3f dyn/cm2 unrelaxed",
      (sgn == +1 ? "lo" : "hi"), qr.p, qf.p);
    check(qr.p < qf.p, "outflow sigma pulls pressure toward target", buf);
  }

  // --- inflow: interior normal velocity above target -> ghost moves down
  //     toward it; interior T above target -> ghost T moves down.
  for (int sgn : {+1, -1}) {
    Real rho = 0.0, e = 0.0;
    eos.PYT2RE(cs.p0, Y, 400.0, rho, e);
    const Real u_int = (sgn == +1) ? 4.0e3 : -4.0e3; // inflow through the face
    Real sN[NVAR], sg[NVAR];
    set_state(sN, rho, u_int, 400.0, Y);
    pc_nscbc::Target t;
    t.type = pc_nscbc::Type::inflow;
    t.u[0] = (sgn == +1) ? 2.0e3 : -2.0e3; // want less inflow
    t.T = 300.0;                           // want colder
    for (int k = 0; k < NUM_SPECIES; k++) {
      t.Y[k] = Y[k];
    }
    pc_nscbc::apply(sN, sN, sN, 1, dx, 0, sgn, 1, t, prm, sg);
    const Prim q = get_prim(sg);
    // "toward the target" in the outward frame
    const Real n_sgn = -static_cast<Real>(sgn);
    const Real uo_int = n_sgn * u_int, uo_g = n_sgn * q.u,
               uo_t = n_sgn * t.u[0];
    char buf[200];
    std::snprintf(
      buf, sizeof(buf), "%s face: u_out %.2f -> %.2f cm/s (target %.2f)",
      (sgn == +1 ? "lo" : "hi"), uo_int, uo_g, uo_t);
    check(
      (uo_g - uo_t) * (uo_int - uo_t) >= 0.0 &&
        std::abs(uo_g - uo_t) < std::abs(uo_int - uo_t),
      "inflow relax_u moves velocity toward target", buf);
    std::snprintf(
      buf, sizeof(buf), "%s face: T %.3f -> %.3f (target %.1f)",
      (sgn == +1 ? "lo" : "hi"), 400.0, q.T, t.T);
    check(
      q.T < 400.0 && q.T > t.T,
      "inflow relax_t moves temperature toward target", buf);
  }
}

// C3: species handling.  At an outflow a composition gradient must be
//     extrapolated, NOT clamped to any target -- the legacy code imposed the
//     target composition at outflows, which over-specifies the problem.  At an
//     inflow the target composition must be imposed exactly.
void
check_species()
{
  Case cs;
  auto eos = pele::physics::PhysicsType::eos();
  const Real dx = cs.L / cs.n;
  pc_nscbc::Params prm;
  prm.L_ref = cs.L;
  prm.sigma = 0.25;

  // A linear composition ramp running out through a hi outflow face.
  Real sN[NVAR], sNm1[NVAR], sNm2[NVAR], sg[NVAR];
  auto mk = [&](Real* s, Real yO2) {
    Real Y[NUM_SPECIES];
    Y[O2_ID] = yO2;
    Y[N2_ID] = 1.0 - yO2;
    Real rho = 0.0, e = 0.0;
    eos.PYT2RE(cs.p0, Y, cs.T0, rho, e);
    set_state(s, rho, 3.0e3, cs.T0, Y);
  };
  mk(sN, 0.30);
  mk(sNm1, 0.28);
  mk(sNm2, 0.26);
  pc_nscbc::Target t;
  t.type = pc_nscbc::Type::outflow;
  t.p = cs.p0;
  pc_nscbc::apply(sN, sNm1, sNm2, 3, dx, 0, -1, 1, t, prm, sg);
  const Prim q = get_prim(sg);
  char buf[200];
  std::snprintf(
    buf, sizeof(buf), "Y(O2) 0.26,0.28,0.30 -> ghost %.6f (expect 0.32)",
    q.Y[O2_ID]);
  check(
    std::abs(q.Y[O2_ID] - 0.32) < 1e-6,
    "outflow extrapolates composition (not imposed)", buf);

  // Sum of mass fractions must be exactly one after extrapolation+clipping.
  Real ysum = 0.0;
  for (int k = 0; k < NUM_SPECIES; k++) {
    ysum += sg[UFS + k];
  }
  std::snprintf(
    buf, sizeof(buf), "|sum(rhoY)/rho - 1| = %.3e",
    std::abs(ysum / sg[URHO] - 1.0));
  check(
    std::abs(ysum / sg[URHO] - 1.0) < 1e-14, "sum(Y) == 1 after outflow fill",
    buf);

  // Inflow imposes the target exactly.
  pc_nscbc::Target ti;
  ti.type = pc_nscbc::Type::inflow;
  ti.T = cs.T0;
  ti.u[0] = -3.0e3;
  ti.Y[O2_ID] = 1.0;
  ti.Y[N2_ID] = 0.0;
  Real sIn[NVAR];
  mk(sIn, 0.30);
  set_state(sIn, get_prim(sIn).rho, -3.0e3, cs.T0, [&] {
    static Real Y[NUM_SPECIES];
    Y[O2_ID] = 0.30;
    Y[N2_ID] = 0.70;
    return Y;
  }());
  pc_nscbc::apply(sIn, sIn, sIn, 1, dx, 0, -1, 1, ti, prm, sg);
  const Prim qi = get_prim(sg);
  std::snprintf(
    buf, sizeof(buf), "interior Y(O2)=0.30, target 1.0 -> ghost %.6f",
    qi.Y[O2_ID]);
  check(
    std::abs(qi.Y[O2_ID] - 1.0) < 1e-12, "inflow imposes composition exactly",
    buf);

  // Energy identity must hold exactly.
  const Real ke = 0.5 * sg[UMX] * sg[UMX] / sg[URHO];
  const Real resid =
    std::abs(sg[UEDEN] - sg[UEINT] - ke) / std::max(std::abs(sg[UEDEN]), 1.0);
  std::snprintf(buf, sizeof(buf), "|UEDEN - UEINT - KE|/UEDEN = %.3e", resid);
  check(resid < 1e-14, "UEDEN == UEINT + KE exactly", buf);
}

// C4: acoustic reflection.  Launch a Gaussian pressure pulse at an outflow and
//     measure the amplitude that comes back.  Reports the sigma sweep, which
//     is the curve the AMReX regression test is checked against.
Real
reflection_coefficient(Real sigma, int n, int order, bool pin, Real* p_drift)
{
  Case cs;
  cs.n = n;
  Field f(cs.n), g(cs.n), h(cs.n);
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
  auto eos = pele::physics::PhysicsType::eos();
  const Real dx = cs.L / cs.n;

  const Real amp = 1.0e-3; // 0.1% pressure pulse -- linear regime
  const Real x0 = 0.35 * cs.L, w = 0.04 * cs.L;
  Real rho_ref = 0.0, e_ref = 0.0;
  eos.PYT2RE(cs.p0, Y, cs.T0, rho_ref, e_ref);
  const Prim qref = [&] {
    Prim q{};
    q.rho = rho_ref;
    q.u = 0.0;
    q.p = cs.p0;
    q.T = cs.T0;
    for (int k = 0; k < NUM_SPECIES; k++) {
      q.Y[k] = Y[k];
    }
    return q;
  }();
  const Real c0 = sound_speed(qref);

  for (int i = -NG; i < cs.n + NG; i++) {
    const Real x = (i + 0.5) * dx;
    const Real dp = amp * cs.p0 * std::exp(-std::pow((x - x0) / w, 2));
    // Right-running isentropic acoustic pulse.
    const Real p = cs.p0 + dp;
    const Real rho = rho_ref + dp / (c0 * c0);
    const Real u = dp / (rho_ref * c0);
    Real T = 0.0;
    eos.RYP2T(rho, Y, p, T);
    set_state(f.at(i), rho, u, T, Y);
  }

  pc_nscbc::Params prm;
  prm.L_ref = cs.L;
  prm.sigma = sigma;
  prm.order = order;
  prm.pin_farfield = pin;
  pc_nscbc::Target off, out;
  out.type = pc_nscbc::Type::outflow;
  out.p = cs.p0;

  const Real t_end = 1.6 * cs.L / c0; // pulse out, any reflection back in
  Real t = 0.0;
  Real peak_in = amp * cs.p0, peak_out = 0.0;
  // measurement window: the left half, after the pulse has left
  const Real t_measure = 0.9 * cs.L / c0;

  while (t < t_end) {
    Real cmax = 0.0;
    for (int i = 0; i < cs.n; i++) {
      const Prim q = get_prim(f.at(i));
      cmax = std::max(cmax, std::abs(q.u) + sound_speed(q));
    }
    Real dt = cs.cfl * dx / cmax;
    dt = std::min(dt, t_end - t);
    if (dt <= 0.0) {
      break;
    }

    fill_bcs(f, off, out, prm, dx);
    stage(f, g, dx, dt);
    fill_bcs(g, off, out, prm, dx);
    stage(g, h, dx, dt);
    for (int i = 0; i < cs.n; i++) {
      for (int v = 0; v < NVAR; v++) {
        f.at(i)[v] = 0.5 * (f.at(i)[v] + h.at(i)[v]);
      }
    }
    t += dt;

    if (t > t_measure) {
      // Measure the wave content only: the deviation from the instantaneous
      // domain mean.  Measuring against p0 instead would count the sigma
      // anchoring transient (a uniform, non-propagating pressure adjustment)
      // as if it were a reflected wave, which it is not.
      Real pbar = 0.0;
      for (int i = 0; i < cs.n; i++) {
        pbar += get_prim(f.at(i)).p;
      }
      pbar /= cs.n;
      for (int i = 0; i < cs.n / 2; i++) {
        peak_out = std::max(peak_out, std::abs(get_prim(f.at(i)).p - pbar));
      }
    }
  }
  if (p_drift != nullptr) {
    Real pm = 0.0;
    for (int i = 0; i < cs.n; i++) {
      pm += get_prim(f.at(i)).p;
    }
    *p_drift = pm / cs.n - cs.p0;
  }
  return peak_out / peak_in;
}

void
check_reflection(bool sweep)
{
  char buf[220];
  Real drift = 0.0;
  const Real R2 = reflection_coefficient(0.25, 400, 2, false, &drift);
  std::snprintf(
    buf, sizeof(buf), "R = %.4f %% (sigma=0.25, n=400, order=2)", 100.0 * R2);
  check(R2 < 0.01, "outflow reflection below 1%", buf);

  // Judge the extrapolation order at sigma = 0.  At sigma = 0.25 the residual
  // in the domain is dominated by the anchoring transient, not by reflection,
  // so an order comparison there measures the wrong thing.
  const Real R0o1 = reflection_coefficient(0.0, 400, 1, false, nullptr);
  const Real R0o2 = reflection_coefficient(0.0, 400, 2, false, nullptr);
  std::snprintf(
    buf, sizeof(buf), "sigma=0: order1 R = %.5f %%, order2 R = %.5f %%",
    100.0 * R0o1, 100.0 * R0o2);
  check(
    R0o2 <= R0o1 * 1.05, "2nd-order extrapolation is not worse than 1st", buf);

  std::snprintf(
    buf, sizeof(buf), "sigma=0 R = %.5f %%, sigma=0.25 R = %.5f %%",
    100.0 * R0o2, 100.0 * R2);
  check(R0o2 <= R2 * 1.05, "sigma=0 is the least reflecting", buf);

  if (sweep) {
    std::printf("\n  sigma sweep (n=400, order=2)\n");
    std::printf(
      "  %8s  %12s  %16s  %14s\n", "sigma", "R [%]", "p drift [dyn/cm2]",
      "tau_relax [s]");
    for (Real s : {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1.0, 2.0}) {
      Real d = 0.0;
      const Real R = reflection_coefficient(s, 400, 2, false, &d);
      const Real tau = (s > 0.0) ? Case().L / (s * 34783.7) : 0.0;
      std::printf("  %8.3f  %12.5f  %16.6e  %14.4e\n", s, 100.0 * R, d, tau);
    }
    Real d = 0.0;
    const Real Rp = reflection_coefficient(0.0, 400, 2, true, &d);
    std::printf(
      "  %8s  %12.5f  %16.6e  %14s   (pin_farfield)\n", "--", 100.0 * Rp, d,
      "0 (value pin)");
  }
}

// C5: the relaxation rate must be a RATE -- grid-independent, and equal to
//     K = sigma (1-M^2) c / L.  This is the check that distinguishes the
//     Poinsot-Lele parameterisation adopted here from CAMR's value-blend,
//     whose effective rate is c/dx and therefore doubles when the mesh does.
void
check_relaxation_rate()
{
  Case cs;
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
  auto eos = pele::physics::PhysicsType::eos();
  const Real sigma = 0.5;

  auto measure = [&](int n) {
    Case c2;
    c2.n = n;
    Field f(n), g(n), h(n);
    const Real dx = c2.L / n;
    const Real dp0 = 0.02 * c2.p0;
    Real rho = 0.0, e = 0.0, T = 0.0;
    eos.PYT2RE(c2.p0 + dp0, Y, c2.T0, rho, e);
    eos.RYP2T(rho, Y, c2.p0 + dp0, T);
    for (int i = -NG; i < n + NG; i++) {
      set_state(f.at(i), rho, 0.0, T, Y);
    }
    const Prim q0 = get_prim(f.at(0));
    const Real c0 = sound_speed(q0);

    pc_nscbc::Params prm;
    prm.L_ref = c2.L;
    prm.sigma = sigma;
    pc_nscbc::Target out;
    out.type = pc_nscbc::Type::outflow;
    out.p = c2.p0;
    // Both faces outflow so the box simply depressurises.
    pc_nscbc::Target out2 = out;

    const Real t_end = 3.0 * c2.L / c0;
    Real t = 0.0;
    while (t < t_end) {
      Real cmax = 0.0;
      for (int i = 0; i < n; i++) {
        const Prim q = get_prim(f.at(i));
        cmax = std::max(cmax, std::abs(q.u) + sound_speed(q));
      }
      Real dt = std::min(c2.cfl * dx / cmax, t_end - t);
      if (dt <= 0.0) {
        break;
      }
      fill_bcs(f, out2, out, prm, dx);
      stage(f, g, dx, dt);
      fill_bcs(g, out2, out, prm, dx);
      stage(g, h, dx, dt);
      for (int i = 0; i < n; i++) {
        for (int v = 0; v < NVAR; v++) {
          f.at(i)[v] = 0.5 * (f.at(i)[v] + h.at(i)[v]);
        }
      }
      t += dt;
    }
    Real pm = 0.0;
    for (int i = 0; i < n; i++) {
      pm += get_prim(f.at(i)).p;
    }
    pm /= n;
    // exp decay over t_end
    const Real ratio = std::max((pm - c2.p0) / dp0, 1e-12);
    return std::make_pair(-std::log(ratio) / t_end, c0);
  };

  const auto r200 = measure(200);
  const auto r800 = measure(800);
  const Real K_expected = sigma * r200.second / Case().L;
  char buf[240];
  std::snprintf(
    buf, sizeof(buf), "K(n=200)=%.2f, K(n=800)=%.2f 1/s  (ratio %.3f)",
    r200.first, r800.first, r800.first / r200.first);
  check(
    std::abs(r800.first / r200.first - 1.0) < 0.15,
    "relaxation rate is grid-independent", buf);
  std::snprintf(
    buf, sizeof(buf), "measured %.2f vs sigma*c/L = %.2f 1/s (ratio %.3f)",
    r800.first, K_expected, r800.first / K_expected);
  check(
    r800.first / K_expected > 0.2 && r800.first / K_expected < 5.0,
    "relaxation rate is within an order of K=sigma*c/L", buf);
}

// C6: robustness.  Every fallback path must return a finite, physical state
//     and must be counted.
void
check_fallbacks()
{
  Case cs;
  Real Y[NUM_SPECIES];
  for (int k = 0; k < NUM_SPECIES; k++) {
    Y[k] = air_Y(k);
  }
  auto eos = pele::physics::PhysicsType::eos();
  const Real dx = cs.L / cs.n;
  pc_nscbc::Params prm;
  prm.L_ref = cs.L;
  pc_nscbc::Target out;
  out.type = pc_nscbc::Type::outflow;
  out.p = cs.p0;

  Real rho = 0.0, e = 0.0;
  eos.PYT2RE(cs.p0, Y, cs.T0, rho, e);

  int diag[pc_nscbc::Diag::count] = {0};
  Real sN[NVAR], sg[NVAR];

  // supersonic outflow
  set_state(sN, rho, 8.0e4, cs.T0, Y);
  pc_nscbc::apply(sN, sN, sN, 1, dx, 0, -1, 1, out, prm, sg, diag);
  bool finite = true;
  for (int v = 0; v < NVAR; v++) {
    finite = finite && std::isfinite(sg[v]);
  }
  check(
    finite && diag[pc_nscbc::Diag::supersonic] == 1 && sg[URHO] == sN[URHO],
    "supersonic outflow -> exact zero-gradient copy", "counted, state finite");

  // reversed flow at an outflow face
  set_state(sN, rho, -3.0e3, cs.T0, Y);
  pc_nscbc::apply(sN, sN, sN, 1, dx, 0, -1, 1, out, prm, sg, diag);
  finite = true;
  for (int v = 0; v < NVAR; v++) {
    finite = finite && std::isfinite(sg[v]);
  }
  check(
    finite && diag[pc_nscbc::Diag::reversed] == 1 && sg[URHO] > 0.0,
    "reversed flow at outflow -> ambient closure", "counted, state finite");

  // EB body state in the stencil
  Real body[NVAR];
  for (int v = 0; v < NVAR; v++) {
    body[v] = -1.0;
  }
  set_state(sN, rho, 3.0e3, cs.T0, Y);
  const int before = diag[pc_nscbc::Diag::body_state];
  pc_nscbc::apply(sN, body, body, 3, dx, 0, -1, 1, out, prm, sg, diag);
  finite = true;
  for (int v = 0; v < NVAR; v++) {
    finite = finite && std::isfinite(sg[v]);
  }
  check(
    finite && diag[pc_nscbc::Diag::body_state] > before && sg[URHO] > 0.0,
    "covered cells in stencil -> order degraded, no FPE",
    "counted, state finite");

  // fully covered boundary cell
  pc_nscbc::apply(body, body, body, 3, dx, 0, -1, 1, out, prm, sg, diag);
  finite = true;
  for (int v = 0; v < NVAR; v++) {
    finite = finite && std::isfinite(sg[v]);
  }
  check(finite, "covered boundary cell -> finite fallback", "state finite");
}

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv, false);
  {
    const bool sweep = (argc > 1) && (std::string(argv[1]) == "sweep");
    std::printf("\nnscbc1d -- standalone verification of Source/NSCBC.H\n");
    std::printf(
      "EOS: %s,  NUM_SPECIES = %d,  NVAR = %d\n\n",
      pele::physics::PhysicsType::identifier().c_str(), NUM_SPECIES, NVAR);

    std::printf("C1  uniform-state consistency\n");
    check_uniform();
    std::printf("\nC2  relaxation directions\n");
    check_relaxation_signs();
    std::printf("\nC3  species and state identities\n");
    check_species();
    std::printf("\nC4  acoustic reflection\n");
    check_reflection(sweep);
    std::printf("\nC5  relaxation rate is a rate\n");
    check_relaxation_rate();
    std::printf("\nC6  fallbacks\n");
    check_fallbacks();

    std::printf("\n%d passed, %d failed\n\n", n_pass, n_fail);
  }
  amrex::Finalize();
  return (n_fail == 0) ? 0 : 1;
}
