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

template<class Tpsys>
void check_val(const Tpsys& sys, const int tag) {
  std::cout << "tag is " << tag << std::endl;
  for(PS::S32 i = 0; i < sys.getNumberOfParticleLocal(); i++) {
    assert(std::isfinite(sys[i].acc.x));
    assert(std::isfinite(sys[i].acc.y));
    assert(std::isfinite(sys[i].acc.z));

    assert(std::isfinite(sys[i].press.x));
    assert(std::isfinite(sys[i].press.y));
    assert(std::isfinite(sys[i].press.z));
  }
}

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

  if (PS::Comm::getRank() == 0) param.CalcCorePtclId(system);

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

  std::cout << "Num / proc is " << system.getNumberOfParticleLocal() << std::endl;;

  PS::TreeForForceShort<RESULT::ForceDPD, EPI::DPD, EPJ::DPD>::Gather force_tree;
  force_tree.initialize(3 * system.getNumberOfParticleGlobal() );
  force_tree.calcForceAllAndWriteBack(CalcForceEpEpDPD(), system, dinfo);

  check_val(system, 0);

  for (PS::S32 i = 0; i < system.getNumberOfParticleLocal(); i++) {
    system[i].acc = 0.0;
  }

  // bonded test
  PS::F64vec bonded_vir(0.0, 0.0, 0.0);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  const PS::U32 est_loc_amp = param.init_amp_num / PS::Comm::getNumberOfProc() + 100;
  ForceBondedMPI<PS::ParticleSystem<FPDPD>, EPJ::DPD> fbonded(est_loc_amp);
  fbonded.CalcListedForce(system, force_tree.epj_org(), bonded_vir, param.core_amp_id());
#else
  ForceBonded<PS::ParticleSystem<FPDPD> > fbonded(system, Parameter::all_unit * param.init_amp_num, param.init_amp_num, param.core_amp_id());
  fbonded.CalcListedForce(system, bonded_vir);
#endif

  check_val(system, 1);
  
  Observer<PS::ParticleSystem<FPDPD> > observer(cdir);
  observer.Initialize();
  observer.Trajectory(system);
  observer.Pressure(system, bonded_vir, Parameter::ibox_leng);
  
  PS::Finalize();
}
