#include <iostream>
#include "particle_simulator.hpp"
#include "io_util.hpp"
#include "parameter.hpp"
#include "f_calculator.hpp"
#include "observer.hpp"

constexpr char Parameter::atom_type[21];

PS::F64vec Parameter::box_leng, Parameter::ibox_leng;
PS::U32 Parameter::time;
PS::U32 Parameter::all_time, Parameter::step_mic, Parameter::step_mac;

int main(int argc, char *argv[]) {
  PS::Initialize(argc, argv);
  if (argc != 2) {
    std::cerr << "argv[1] = target directory.\n";
    PS::Abort();
  }

  const std::string cdir = argv[1];
  Parameter param(cdir);
  param.Initialize();
  if (PS::Comm::getRank() == 0) param.LoadParam();

  PS::ParticleSystem<FPDPD> system;
  system.initialize();
  if(PS::Comm::getRank() == 0) {
    system.setNumberOfParticleLocal(param.init_amp_num * Parameter::all_unit + param.sol_num);
    Parameter::time = param.LoadParticleConfig(system);
  } else {
    system.setNumberOfParticleLocal(0);
  }

  if (PS::Comm::getRank() == 0) param.CalcCorePtclId(system);

  param.ShareDataWithOtherProc();
  param.CheckLoaded();
  param.CheckParticleConfigIsValid(system);
  param.DebugDumpAllParam();

  PS::DomainInfo dinfo;
  const PS::F64 coef_ema = 0.3;
  dinfo.initialize(coef_ema);
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), Parameter::box_leng);
  dinfo.decomposeDomainAll(system);
  system.exchangeParticle(dinfo);

  std::cout << system.getNumberOfParticleLocal() << std::endl;;

  Observer<PS::ParticleSystem<FPDPD> > observer(cdir);
  observer.Initialize();
  observer.Trajectory(system);

  PS::Finalize();
}
