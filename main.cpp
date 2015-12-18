#include <iostream>
#include <fstream>
#include "f_calculator.hpp"
#include "driftkick.hpp"
#include "parameter.hpp"
#include "observer.hpp"
#include <sys/time.h>

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#error "MPI is not supported yet!"
#endif

namespace  {
  timeval tv_beg, tv_end;
  double timme_diff(timeval& tv1, timeval& tv2) {
    return tv2.sec - tv1.sec + (tv2.usec - tv1.usec) * 1.0e-6;
  }
  void timer_start() { gettimeofday(&tv_beg, NULL); }
  void timer_stop () { gettimeofday(&tv_end, NULL); }
  void show_duration() {
    int all_sec = static_cast<int>(time_diag(tv_beg, tv_end));
    const int d = all_sec / (60 * 60 * 24);
    all_sec %= 60 * 60 * 24;
    const int h = all_sec / (60 * 60);
    all_sec %= 60 * 60;
    const int m = all_sec / 60;
    all_sec %= 60;
    const int s = all_sec;
    std::cout << "Simulation time is \n";
    std::cout << d << "d " << h << "h " << m << "m " << "s \n";
  }
}

int main(int argc, char *argv[]) {
  timer_start();
  
  PS::Initialize(argc, argv);
  
  Parameter param(std::string(argv[1]));
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
  dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), param.box_leng);
  dinfo.collectSampleParticle(system);
  dinfo.decomposeDomain();
  system.exchangeParticle(dinfo);
  
  PS::TreeForForceShort<ForceDPD, EPIDPD, EPJDPD>::Scatter tree;
  tree.initialize(3 * system.getNumberOfParticleGlobal() );
  tree.calcForceAllAndWriteBack(CalcForceEpEp(), system, dinfo);

  //observer
  Observer<PS::ParticleSystem<FPDPD> > observer;
  observer.Initialize();

  //main loop
  const int all_time = 100000, step_mic = 1000, step_mac = 100;
  for(Parameter::time = 0; time < all_time; time++) {
    kick_and_drift(system, param.dt, param.box_leng, param.ibox_leng);
    
    tree.initialize(3 * system.getNumberOfParticleGlobal()); //NOTE: this routine may not be called.
    dinfo.decomposeDomain();
    system.exchangeParticle(dinfo);
    tree.calcForceAllAndWriteBack(CalcForceEpEp(), system, dinfo);
    
    kick(system, param.dt);

    //observer
    if(time % step_mac == 0) {
      observer.KineticTempera(system);
      observer.Pressure(system, param.ibox_leng);
      //observer.ConfigTempera();
      //observer.Diffusion();
      
    }

    if(time % time_mic == 0) {
      observer.Trajectory(system);
      
    }
    
  }//end of main loop
  
  timer_stop();
  show_duration();
  
  PS::Finalize();
}
