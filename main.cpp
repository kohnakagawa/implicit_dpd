#include <iostream>
#include <fstream>
#include "f_calculator.hpp"
#include "driftkick.hpp"
#include "parameter.hpp"

int main(int argc, char *argv[]) {
  PS::F64 Tbegin = PS::GetWtime();
  std::cout << std::setprecision(15);
  std::cerr << std::setprecision(15);
  
  PS::Initialize(argc, argv);
  
  Parameter param(std::string(argv[1]));
  param.Initialize();
  param.LoadParam();
  param.LoadCheck();
  
  //copy interactions
  for(int i = 0; i < Parameter::prop_num; i++) {
    for(int j = 0; j < Parameter::prop_num; j++) {
      CalcForceEpEp::cf_c[i][j] = param.cf_c[i][j];
      CalcForceEpEp::cf_g[i][j] = param.cf_g[i][j];
      CalcForceEpEp::cf_r[i][j] = param.cf_r[i][j];
    }
  }


  PS::ParticleSystem<FPDPD> system;
  system.initialize();
  //load particle configuration.

  //or make particle configuration.
  
  
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

  //main loop
  const int all_time = 100000, step_mic = 1000, step_mac = 100;
  for(int time = 0; time < all_time; time++) {
    KickAndDrift(system, param.dt, param.box_leng, param.ibox_leng);
    
    tree.initialize(3 * system.getNumberOfParticleGlobal()); //NOTE this routine may not be called.
    dinfo.decomposeDomain();
    system.exchangeParticle(dinfo);
    tree.calcForceAllAndWriteBack(CalcForceEpEp(), system, dinfo);
    
    Kick(system, param.dt);

    //observer
    
  }//end of main loop
  
  PS::Finalize();
}
