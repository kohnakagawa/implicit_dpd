#include <iostream>
#include "particle_simulator.hpp"
#include "io_util.hpp"
#include "parameter.hpp"
#include "f_calc.hpp"
#include "observer.hpp"
#include "driftkick.hpp"
#include <sys/time.h>

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#error "MPI is not supported yet!"
#endif

static_assert(Parameter::head_unit == 1, "head_unit should be 1.");
static_assert(Parameter::tail_unit == 3, "tail_unit should be 3.");
static_assert(Parameter::bond_leng != 0.0, "bond_leng should be not 0.0.");
static_assert(Parameter::Reo < 3.0, "Reo should be less than 3.0.");

namespace  {
  timeval tv_beg, tv_end;
  double time_diff(timeval& tv1, timeval& tv2) {
    return tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec) * 1.0e-6;
  }
  void timer_start() { gettimeofday(&tv_beg, NULL); }
  void timer_stop () { gettimeofday(&tv_end, NULL); }
  void show_duration() {
    int all_sec = static_cast<int>(time_diff(tv_beg, tv_end));
    const int d = all_sec / (60 * 60 * 24);
    all_sec %= 60 * 60 * 24;
    const int h = all_sec / (60 * 60);
    all_sec %= 60 * 60;
    const int m = all_sec / 60;
    all_sec %= 60;
    const int s = all_sec;
    std::cout << "Simulation time is \n";
    std::cout << d << "d " << h << "h " << m << "m " << s << "s\n";
  }
  
  template<class Tpsys>
  void set_local_id(Tpsys& sys) {
    const PS::S32 num = sys.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < num; i++)
      sys[i].id_loc = i;
  }
}

int main(int argc, char *argv[]) {

  
  PS::Initialize(argc, argv);
  
  const std::string cdir = argv[1];
  Parameter param(cdir);
  param.Initialize();
  param.LoadParam();
  param.CheckLoaded();

  PS::ParticleSystem<FPDPD> system;
  system.initialize();
  system.setNumberOfParticleLocal(param.init_amp_num * Parameter::all_unit);
  //load particle configuration.
  param.LoadParticleConfig(system);
  param.CheckParticleConfigIsValid(system);

  PS::DomainInfo dinfo;
  const PS::F64 coef_ema = 0.3;
  dinfo.initialize(coef_ema);
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), Parameter::box_leng);
  dinfo.collectSampleParticle(system);
  dinfo.decomposeDomain();
  system.exchangeParticle(dinfo);

#ifdef USE_NEIGH_BUFFER
  NearlistBuffer n_buffer;
  n_buffer.Initialize(3 * system.getNumberOfParticleLocal() );
  set_local_id(system);

  PS::TreeForForceShort<RESULT::ForceDPDDensity, EPI::DPD_Density, EPJ::DPD_Density>::Gather force_tree;
  force_tree.initialize(3 * system.getNumberOfParticleGlobal() );
  force_tree.calcForceAllAndWriteBack(CalcDensity(), system, dinfo);
  force_tree.calcForceCopyPtclAndWriteBack(CalcForceEpEpDPD(), system, dinfo);

#else
  PS::TreeForForceShort<RESULT::Density, EPI::Density, EPJ::Density>::Gather dens_tree;
  dens_tree.initialize(3 * system.getNumberOfParticleGlobal() );
  dens_tree.calcForceAllAndWriteBack(CalcDensity(), system, dinfo);
  
  PS::TreeForForceShort<RESULT::ForceDPD, EPI::DPD, EPJ::DPD>::Gather force_tree;
  force_tree.initialize(3 * system.getNumberOfParticleGlobal() );
  force_tree.calcForceAllAndWriteBack(CalcForceEpEpDPD(), system, dinfo);
#endif
  
  ForceBonded<PS::ParticleSystem<FPDPD> > fbonded(system, Parameter::all_unit * param.init_amp_num);
  fbonded.CalcListedForce(system);

  // observer
  Observer<PS::ParticleSystem<FPDPD> > observer(cdir);
  observer.Initialize();

  timer_start();  
  //main loop
  const PS::U32 atime = Parameter::all_time;
  PS::F64vec bonded_vir(0.0, 0.0, 0.0);
  for(Parameter::time = 0; Parameter::time < atime; Parameter::time++) {
#ifdef DEBUG    
    std::cout << Parameter::time << std::endl;
#endif
    drift_and_predict(system, param.dt, param.box_leng, param.ibox_leng);
    
    dinfo.decomposeDomain();
    system.exchangeParticle(dinfo);

#ifdef USE_NEIGH_BUFFER
    n_buffer.CheckBufferSize(system.getNumberOfParticleLocal());
    n_buffer.ClearList();
    set_local_id(system);
    //calc density
    force_tree.calcForceAllAndWriteBack(CalcDensity(), system, dinfo);
    
    //calc non-bonded interaction
    //NOTE: After call this function, FPDPD.density information is set to 0.0 .
    force_tree.calcForceCopyPtclAndWriteBack(CalcForceEpEpDPD(), system, dinfo);
#else
    dens_tree.calcForceAllAndWriteBack(CalcDensity(), system, dinfo);
    force_tree.calcForceAllAndWriteBack(CalcForceEpEpDPD(), system, dinfo);
#endif

    fbonded.CalcListedForce(system, bonded_vir);
    
    kick(system, param.dt);

    if(Parameter::time % Parameter::step_mac == 0) {
      observer.KineticTempera(system);
      observer.Pressure(system, bonded_vir, param.ibox_leng);
      observer.Diffusion(system, param.amp_num);
      //observer.ConfigTempera();
    }

    if(Parameter::time % Parameter::step_mic == 0) {
      observer.Trajectory(system);
    }

#ifdef DEBUG
    if(Parameter::time % Observer<PS::ParticleSystem<FPDPD> >::flush_freq == 0)
      observer.FlushAll();
#endif
    
  }//end of main loop
  
  timer_stop();
  show_duration();

  //print configuration for restart
  observer.FinConfig(system);
  
  observer.CleanUp();
  param.DumpAllParam();
  
  PS::Finalize();
}