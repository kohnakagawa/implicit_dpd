#include <iostream>
#include "particle_simulator.hpp"
#include "parameter.hpp"
#include "f_calculator.hpp"
#include "observer.hpp"
#include "driftkick.hpp"
#include <sys/time.h>

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#error "MPI is not supported yet!"
#endif

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
}

int main(int argc, char *argv[]) {
  timer_start();
  
  PS::Initialize(argc, argv);
  
  const std::string cdir = argv[1];
  Parameter param(cdir);
  param.Initialize();
  param.LoadParam();
  param.CheckLoaded();
  
  PS::ParticleSystem<FPDPD> system;
  system.initialize();
  system.setNumberOfParticleLocal(param.init_prtcl_num);
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
  
  PS::TreeForForceShort<ForceDPD, EPIDPD, EPJDPD>::Scatter tree;
  tree.initialize(3 * system.getNumberOfParticleGlobal() );
  tree.calcForceAllAndWriteBack(CalcForceEpEp(), system, dinfo);

  //observer
  Observer<PS::ParticleSystem<FPDPD> > observer(cdir);
  observer.Initialize();

  //main loop
  const PS::U32 all_time = 100000, step_mic = 1000, step_mac = 100;
  for(Parameter::time = 0; Parameter::time < all_time; Parameter::time++) {
    kick_and_drift(system, param.dt, param.box_leng, param.ibox_leng);
    
    tree.initialize(3 * system.getNumberOfParticleGlobal()); //NOTE: this routine may not be called.
    dinfo.decomposeDomain();
    system.exchangeParticle(dinfo);
    tree.calcForceAllAndWriteBack(CalcForceEpEp(), system, dinfo);
    
    kick(system, param.dt);

    //observer
    if(Parameter::time % step_mac == 0) {
      observer.KineticTempera(system);
      observer.Pressure(system, param.ibox_leng);
      //observer.ConfigTempera();
      //observer.Diffusion();
      
    }

    if(Parameter::time % step_mic == 0) {
      observer.Trajectory(system);
    }
    
  }//end of main loop
  
  timer_stop();
  show_duration();
  
  PS::Finalize();
}
