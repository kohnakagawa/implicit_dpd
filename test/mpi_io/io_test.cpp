#include <iostream>
#include "particle_simulator.hpp"
#include "io_util.hpp"
#include "parameter.hpp"
#include "f_calculator.hpp"
#include "observer.hpp"

int main(int argc, char *argv[]) {
  PS::Initialize(argc, argv);

  const std::string cdir = argv[1];  
  Parameter param(cdir);
  param.Initialize();
  if(PS::Comm::getRank() == 0) param.LoadParam();
  
  PS::ParticleSystem<FPDPD> system;
  system.initialize();
  if(PS::Comm::getRank() == 0) {
    system.setNumberOfParticleLocal(param.init_amp_num * Parameter::all_unit);
    Parameter::time = param.LoadParticleConfig(system);
  } else {
    system.setNumberOfParticleLocal(0);    
  }

  param.ShareDataWithOtherProc();
  param.CheckLoaded();
  param.DebugDumpAllParam();

  PS::DomainInfo dinfo;
  const PS::F64 coef_ema = 0.3;
  dinfo.initialize(coef_ema);
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), Parameter::box_leng);
  dinfo.collectSampleParticle(system);
  dinfo.decomposeDomain();
  system.exchangeParticle(dinfo);

  std::cout << system.getNumberOfParticleLocal() << std::endl;;
  
  Observer<PS::ParticleSystem<FPDPD> > observer(cdir);
  observer.Initialize();
  observer.Trajectory(system);
  
  PS::Finalize();
}
