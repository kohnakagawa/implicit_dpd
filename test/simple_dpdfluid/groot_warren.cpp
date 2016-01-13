#include <iostream>
#include "particle_simulator.hpp"
#include "io_util.hpp"
#include "parameter.hpp"
#include "f_calculator.hpp"
#include "observer.hpp"
#include "driftkick.hpp"

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#error "MPI is not supported yet!"
#endif

static_assert(Parameter::all_unit == 4, "all_unit should be 4.");

namespace {
  void calc_mean_pressure(const std::string& cdir,
			  const PS::S32 beg,
			  const PS::F64 rho)
  {
    std::string fname = cdir + "/pressure.txt";
    std::ifstream fin(fname.c_str());
    assert(fin);

    PS::F64vec press_mean(0.0, 0.0, 0.0), buf(0.0, 0.0, 0.0);
    PS::S32 cnt = 0, time = 0;
    
    while(true) {
      fin >> buf;
      if(fin.eof()) break;
      if(time >= beg) {
	press_mean += buf;
	cnt++;
      }
      time++;
    }
    
    press_mean /= cnt;
    const PS::F64 press = (press_mean.x + press_mean.y + press_mean.z) / 3.0;
    fname = "./rho_vs_press.txt";
    std::ofstream fout(fname.c_str(), std::ios::app);
    fout << std::setprecision(15);
    fout << rho << " " << press << " " << (press - rho) / (rho * rho * 25.0) << std::endl; //NOTE: calculate only potential component.
  }

  void calc_mean_kintemp(const std::string& cdir,
			 const PS::S32 beg,
			 const PS::F64 dt)
  {
    std::string fname = cdir + "/kin_temp.txt";
    std::ifstream fin(fname.c_str());
    assert(fin);

    PS::F64 T_mean = 0.0, buf = 0.0;
    PS::S32 cnt = 0, time = 0;
    
    while(true) {
      fin >> buf;
      if(fin.eof()) break;
      if(time >= beg) {
	T_mean += buf;
	cnt++;
      }
      time++;
    }

    T_mean /= cnt;
    fname = "./dt_vs_tempera.txt";
    std::ofstream fout(fname.c_str(), std::ios::app);
    fout << std::setprecision(15);
    fout << dt << " " << T_mean << std::endl;
  }

  enum {
    DRIFT_KICK = 0,
    DECOMPOSE,
    EXCHANGE,
    F_CALC,
    KICK,
  };
  
  template<class Tpsys>
  bool check_value_is_finite(Tpsys& system, const PS::U32 id) {
    const PS::U32 num = system.getNumberOfParticleLocal();
    for(PS::U32 i = 0; i < num; i++) {
      const bool pos_valid = (std::isfinite(system[i].pos.x) &&
			      std::isfinite(system[i].pos.y) &&
			      std::isfinite(system[i].pos.z));
      const bool vel_valid = (std::isfinite(system[i].vel.x) &&
			      std::isfinite(system[i].vel.y) &&
			      std::isfinite(system[i].vel.z));
      const bool acc_valid = (std::isfinite(system[i].acc.x) &&
			      std::isfinite(system[i].acc.y) &&
			      std::isfinite(system[i].acc.z));
      if(!pos_valid)
	std::cerr << "id " << id << " position is not finite value.\n";

      if(!vel_valid)
	std::cerr << "id " << id << " velocity is not finite value.\n";

      if(!acc_valid)
	std::cerr << "id " << id << " force is not finite value.\n";
	
      if( (!pos_valid) || (!vel_valid) || (!acc_valid) ) {
	std::cerr << "id: " << system[i].id << std::endl;
	std::cerr << "pos: " << system[i].pos << std::endl;
	std::cerr << "vel: " << system[i].vel << std::endl;
	std::cerr << "acc: " << system[i].acc << std::endl;
	return false;
      }
    }
    return true;
  }
};

#define CHECK_VAL(obj, id)						\
  if(!check_value_is_finite((obj), (id))) {				\
    std::cerr << "Error occurs at " << __FILE__ << " " << __LINE__ << std::endl; \
    PS::Abort();							\
  }

int main(int argc, char *argv[]) {
  PS::Initialize(argc, argv);
  if(argc != 4) {
    std::cerr << "argv[1] = target directory.\n";
    std::cerr << "argv[2] = observe beg time.\n";
    std::cerr << "argv[3] = run mode.\n";
    PS::Abort();
  }
  
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

  PS::TreeForForceShort<RESULT::Density, EPI::Density, EPJ::Density>::Gather dens_tree;
  dens_tree.initialize(3 * system.getNumberOfParticleGlobal() );
  dens_tree.calcForceAllAndWriteBack(CalcDensity(), system, dinfo);

  PS::TreeForForceShort<RESULT::ForceDPD, EPI::DPD, EPJ::DPD>::Gather force_tree;
  force_tree.initialize(3 * system.getNumberOfParticleGlobal() );
  force_tree.calcForceAllAndWriteBack(CalcForceEpEpDPD(), system, dinfo);
  CHECK_VAL(system, F_CALC);

  // observer
  Observer<PS::ParticleSystem<FPDPD> > observer(cdir);
  observer.Initialize();
  
  //main loop
  PS::F64vec buf(0.0, 0.0, 0.0);
  for(Parameter::time = 0; Parameter::time < Parameter::all_time; Parameter::time++) {
    drift_and_predict(system, param.dt, param.box_leng, param.ibox_leng);
    CHECK_VAL(system, DRIFT_KICK);
    
    dinfo.decomposeDomain();
    CHECK_VAL(system, DECOMPOSE);

    system.exchangeParticle(dinfo);
    CHECK_VAL(system, EXCHANGE);

    dens_tree.calcForceAllAndWriteBack(CalcDensity(), system, dinfo);
    force_tree.calcForceAllAndWriteBack(CalcForceEpEpDPD(), system, dinfo);
    CHECK_VAL(system, F_CALC);
    
    kick(system, param.dt);
    CHECK_VAL(system, KICK);

    if(Parameter::time % Parameter::step_mac == 0) {
      observer.KineticTempera(system);
      observer.Pressure(system, buf, param.ibox_leng);
      observer.Diffusion(system, param.amp_num);
      //observer.ConfigTempera();
    }

    if(Parameter::time % Parameter::step_mic == 0) {
      observer.Trajectory(system);
    }

    if(Parameter::time % Observer<PS::ParticleSystem<FPDPD> >::flush_freq == 0)
      observer.FlushAll();
  }//end of main loop

  observer.FinConfig(system);
  observer.CleanUp();
  param.DumpAllParam();
  
  const PS::S32 beg_time = std::atoi(argv[2]);
  const std::string mode = argv[3];
  
  if(mode == "pressure") {
    const PS::F64 rho = system.getNumberOfParticleGlobal() / (param.box_leng.x * param.box_leng.y * param.box_leng.z);
    calc_mean_pressure(cdir, beg_time, rho);
  } else if(mode == "temperature") {
    calc_mean_kintemp(cdir, beg_time, param.dt);
  } else {
    std::cerr << "Unknown mode.\n";
    PS::Abort();
  }
  
  PS::Finalize();
}
