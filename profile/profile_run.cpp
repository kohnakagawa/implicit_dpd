#include <iostream>
#include "particle_simulator.hpp"
#include "io_util.hpp"
#include "parameter.hpp"
#include "f_calculator.hpp"
#include "observer.hpp"
#include "driftkick.hpp"
#include <numeric>
#include <sys/time.h>

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#error "MPI is not supported yet!"
#endif

static_assert(Parameter::bond_leng == 0.0, "bond_leng should be 0.0.");
static_assert(Parameter::Reo == 3.5, "Reo should be 3.5.");
static_assert(Parameter::arc == 0.9, "arc should be 0.9.");

class Profile {
  FILE* fp = nullptr;
  timeval beg = {0.0, 0.0}, end = {0.0, 0.0};
  std::vector<double> dur_t;
    
public:
  enum {
    DRIFT_KICK = 0,
    DECOMPOSE,
    EXCHANGE,
    F_CALC_NB,
    F_CALC_BO,
    KICK,
    OBSERVE,
      
    FILE_NUM,
  };

  explicit Profile(const char* cur_dir) {
    const std::string fname = std::string(cur_dir) + "/elapse_time.txt";
    fp = io_util::xfopen(fname.c_str(), "w");
    dur_t.resize(FILE_NUM, 0.0);
  }
    
  ~Profile() {
    fclose(fp);
  }
    
  inline void TimerStart() {
    gettimeofday(&beg, nullptr);
  }
    
  inline void TimerStop(const int kind) {
    gettimeofday(&end, nullptr);
    dur_t[kind] += (end.tv_sec - beg.tv_sec) + (end.tv_usec - beg.tv_usec) * 1.0e-6;
  }
    
#define STR(val) #val
    
  void ShowResult() const {
    const double sum_dur = std::accumulate(dur_t.begin(), dur_t.end() + FILE_NUM, 0.0);
    fprintf(fp, "%-20s: %10.4f[s] %10.4f[%%]\n", STR(DRIFT_KICK), dur_t[DRIFT_KICK], dur_t[DRIFT_KICK] / sum_dur * 100.0);
    fprintf(fp, "%-20s: %10.4f[s] %10.4f[%%]\n", STR(DECOMPOSE), dur_t[DECOMPOSE], dur_t[DECOMPOSE] / sum_dur * 100.0);
    fprintf(fp, "%-20s: %10.4f[s] %10.4f[%%]\n", STR(EXCHANGE), dur_t[EXCHANGE], dur_t[EXCHANGE] / sum_dur * 100.0);
    fprintf(fp, "%-20s: %10.4f[s] %10.4f[%%]\n", STR(F_CALC_NB), dur_t[F_CALC_NB], dur_t[F_CALC_NB] / sum_dur * 100.0);
    fprintf(fp, "%-20s: %10.4f[s] %10.4f[%%]\n", STR(F_CALC_BO), dur_t[F_CALC_BO], dur_t[F_CALC_BO] / sum_dur * 100.0);
    fprintf(fp, "%-20s: %10.4f[s] %10.4f[%%]\n", STR(KICK), dur_t[KICK], dur_t[KICK] / sum_dur * 100.0);
    fprintf(fp, "%-20s: %10.4f[s] %10.4f[%%]\n", STR(OBSERVE), dur_t[OBSERVE], dur_t[OBSERVE] / sum_dur * 100.0);
    fprintf(fp, "all_duration: %10.4f[s]", sum_dur);
  }

#undef STR
    
};

int main(int argc, char *argv[]) {
  PS::Initialize(argc, argv);
  if(argc != 3) {
    std::cerr << "argv[1] = target directory.\n";
    std::cerr << "argv[2] = observe beg time.\n";
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

  PS::F64vec bonded_vir(0.0, 0.0, 0.0);
  ForceBonded<PS::ParticleSystem<FPDPD> > fbonded(system, Parameter::all_unit * param.init_amp_num);
  fbonded.CalcListedForce(system, bonded_vir);

  // observer
  Observer<PS::ParticleSystem<FPDPD> > observer(cdir);
  observer.Initialize();

  //for profile
  Profile profile(cdir.c_str());
  
  //main loop
  for(Parameter::time = 0; Parameter::time < Parameter::all_time; Parameter::time++) {

    profile.TimerStart();
    drift_and_predict(system, param.dt, param.box_leng, param.ibox_leng);
    profile.TimerStop(Profile::DRIFT_KICK);

    profile.TimerStart();
    dinfo.decomposeDomain();
    profile.TimerStop(Profile::DECOMPOSE);

    profile.TimerStart();
    system.exchangeParticle(dinfo);
    profile.TimerStop(Profile::EXCHANGE);

    profile.TimerStart();
    dens_tree.calcForceAllAndWriteBack(CalcDensity(), system, dinfo);
    force_tree.calcForceAllAndWriteBack(CalcForceEpEpDPD(), system, dinfo);
    profile.TimerStop(Profile::F_CALC_NB);
    
    profile.TimerStart();
    fbonded.CalcListedForce(system, bonded_vir);
    profile.TimerStop(Profile::F_CALC_BO);
    
    profile.TimerStart();
    kick(system, param.dt);
    profile.TimerStop(Profile::KICK);

    profile.TimerStart();
    if(Parameter::time % Parameter::step_mac == 0) {
      observer.KineticTempera(system);
      observer.Pressure(system, bonded_vir, param.ibox_leng);
      //observer.ConfigTempera();
      //observer.Diffusion();
    }

    if(Parameter::time % Parameter::step_mic == 0) {
      observer.Trajectory(system);
    }

    if(Parameter::time % Observer<PS::ParticleSystem<FPDPD> >::flush_freq == 0)
      observer.FlushAll();
    profile.TimerStop(Profile::OBSERVE);
  }//end of main loop

  //print configuration for restart
  observer.FinConfig(system);
  
  observer.CleanUp();
  param.DumpAllParam();

  profile.ShowResult();
  
  PS::Finalize();
}
