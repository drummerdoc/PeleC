#include "prob.H"

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex::Real* /*problo*/,
  const amrex::Real* /*probhi*/)
{
  {
    amrex::ParmParse pp("prob");
    pp.query("p_amb", PeleC::h_prob_parm_device->p_amb);
    pp.query("T_amb", PeleC::h_prob_parm_device->T_amb);
    pp.query("amp", PeleC::h_prob_parm_device->amp);
    pp.query("x0", PeleC::h_prob_parm_device->x0);
    pp.query("width", PeleC::h_prob_parm_device->width);
  }

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[0] = 1.0;
  amrex::Real rho = 0.0, eint = 0.0, cs = 0.0;
  eos.PYT2RE(
    PeleC::h_prob_parm_device->p_amb, massfrac,
    PeleC::h_prob_parm_device->T_amb, rho, eint);
  eos.RTY2Cs(rho, PeleC::h_prob_parm_device->T_amb, massfrac, cs);
  PeleC::h_prob_parm_device->rho_amb = rho;
  PeleC::h_prob_parm_device->cs_amb = cs;

  amrex::Print() << "  NSCBC-Acoustic: rho_amb = " << rho
                 << " g/cc, c_amb = " << cs << " cm/s\n";
}
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
