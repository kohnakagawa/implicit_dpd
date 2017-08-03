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

constexpr char Parameter::atom_type[21];

PS::F64vec Parameter::box_leng, Parameter::ibox_leng;
PS::U32 Parameter::time;
PS::U32 Parameter::all_time, Parameter::step_mic, Parameter::step_mac;

//NOTE: We check the results from J. Chem. Phys. 132, 155104 (2010) in this code.
static_assert(Parameter::all_unit == 16, "all_unit should be 16.");
static_assert(Parameter::head_unit == 4, "head_unit should be 4.");
static_assert(Parameter::tail_unit == 12, "tail_unit should be 12.");
static_assert(Parameter::bond_leng == 0.0, "bond_leng should be 0.0.");
static_assert(Parameter::Reo == 3.5, "Reo should be 3.5.");
static_assert(Parameter::arc == 0.9, "arc should be 0.9.");

namespace {
  void calc_mean_pressure(const std::string& cdir,
                          const PS::S32 beg,
                          const PS::F64vec& box,
                          const PS::F64 liparea)
  {
    std::string fname = cdir + "/pressure.txt";
    std::ifstream fin(fname.c_str());
    assert(fin);

    PS::F64vec press_mean(0.0, 0.0, 0.0), buf(0.0, 0.0, 0.0);
    PS::F64 stens_mean = 0.0;
    PS::S32 cnt = 0, time = 0;

    while (true) {
      fin >> buf;
      if (fin.eof()) break;
      if (time >= beg) {
        press_mean += buf;
        stens_mean += (buf.y - (buf.x + buf.z) * 0.5) * box.y;
        cnt++;
      }
      time++;
    }

    press_mean /= cnt;
    stens_mean /= cnt;
    fname = "./stens.txt";
    std::ofstream fout(fname.c_str(), std::ios::app);
    fout << std::setprecision(15);
    fout << liparea << " " << stens_mean << " " << press_mean << std::endl;
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

    while (true) {
      fin >> buf;
      if (fin.eof()) break;
      if (time >= beg) {
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
    F_CALC_NB,
    F_CALC_BO,
    KICK,
  };

  template<class Tpsys>
  bool check_value_is_finite(Tpsys& system, const PS::U32 id) {
    const PS::U32 num = system.getNumberOfParticleLocal();
    for (PS::U32 i = 0; i < num; i++) {
      const bool pos_valid = (std::isfinite(system[i].pos.x) &&
                              std::isfinite(system[i].pos.y) &&
                              std::isfinite(system[i].pos.z));
      const bool vel_valid = (std::isfinite(system[i].vel.x) &&
                              std::isfinite(system[i].vel.y) &&
                              std::isfinite(system[i].vel.z));
      const bool acc_valid = (std::isfinite(system[i].acc.x) &&
                              std::isfinite(system[i].acc.y) &&
                              std::isfinite(system[i].acc.z));
      if (!pos_valid)
        std::cerr << "id " << id << " position is not finite value.\n";

      if (!vel_valid)
        std::cerr << "id " << id << " velocity is not finite value.\n";

      if (!acc_valid)
        std::cerr << "id " << id << " force is not finite value.\n";

      if ( (!pos_valid) || (!vel_valid) || (!acc_valid) ) {
        std::cerr << "id: " << system[i].id << std::endl;
        std::cerr << "pos: " << system[i].pos << std::endl;
        std::cerr << "vel: " << system[i].vel << std::endl;
        std::cerr << "acc: " << system[i].acc << std::endl;
        return false;
      }
    }
    return true;
  }

  void change_hyphil_to_hyphob() {
    const double cf_tt = Parameter::cf_c[Parameter::Hyphob][Parameter::Hyphob];
    const double cf_ttt = Parameter::cf_m[Parameter::Hyphob][Parameter::Hyphob][Parameter::Hyphob];
    for (PS::S32 i = 0; i < Parameter::prop_num; i++)
      for (PS::S32 j = 0; j < Parameter::prop_num; j++) {
        Parameter::cf_c[i][j] = cf_tt;
      }
    for (PS::S32 i = 0; i < Parameter::prop_num; i++)
      for (PS::S32 j = 0; j < Parameter::prop_num; j++)
        for (PS::S32 k = 0; k < Parameter::prop_num; k++) {
          Parameter::cf_m[i][j][k] = cf_ttt;
        }
  }
}

#if 0
#define CHECK_VAL(obj, id)                                              \
  if(!check_value_is_finite((obj), (id))) {                             \
    std::cerr << "Error occurs at " << __FILE__ << " " << __LINE__ << std::endl; \
    PS::Abort();                                                        \
  }
#else
#define CHECK_VAL(obj, id)
#endif

int main(int argc, char *argv[]) {
  PS::Initialize(argc, argv);
  if (argc != 3) {
    std::cerr << "argv[1] = target directory.\n";
    std::cerr << "argv[2] = observe beg time.\n";
    PS::Abort();
  }

  const std::string cdir = argv[1];
  Parameter param(cdir);
  param.Initialize();
  if (PS::Comm::getRank() == 0) param.LoadParam();

  //load particle configuration.
  PS::ParticleSystem<FPDPD> system;
  system.initialize();
  if (PS::Comm::getRank() == 0) {
    system.setNumberOfParticleLocal(param.init_amp_num * Parameter::all_unit);
    Parameter::time = param.LoadParticleConfig(system);
  } else {
    system.setNumberOfParticleLocal(0);
  }

  if (PS::Comm::getRank() == 0) param.CalcCorePtclId(system);

  param.ShareDataWithOtherProc();
  param.CheckLoaded();
  param.CheckParticleConfigIsValid(system);

  // change prop hyphil -> hyphob solvent -> hyphob
  change_hyphil_to_hyphob();

  PS::DomainInfo dinfo;
  const PS::F64 coef_ema = 0.3;
  dinfo.initialize(coef_ema);
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), Parameter::box_leng);
  dinfo.decomposeDomainAll(system);
  system.exchangeParticle(dinfo);

  // Initial step
  drift_and_predict(system, param.dt, param.box_leng, param.ibox_leng);
  dinfo.decomposeDomainAll(system);
  system.exchangeParticle(dinfo);
  PS::TreeForForceShort<RESULT::Density, EPI::Density, EPJ::Density>::Gather dens_tree;
  PS::TreeForForceShort<RESULT::ForceDPD, EPI::DPD, EPJ::DPD>::Gather force_tree;
  dens_tree.initialize(3 * system.getNumberOfParticleGlobal() );
  force_tree.initialize(3 * system.getNumberOfParticleGlobal() );

  CalcForceEpEpDPD::m_seed = static_cast<PS::U32>(time(nullptr));
  dens_tree.calcForceAllAndWriteBack(CalcDensity(), system, dinfo);
  force_tree.calcForceAllAndWriteBack(CalcForceEpEpDPD(), system, dinfo);
  CHECK_VAL(system, F_CALC_NB);

  // observer
  Observer<PS::ParticleSystem<FPDPD> > observer(cdir);
  observer.Initialize();
  PS::F64vec bonded_vir(0.0, 0.0, 0.0);

  kick(system, param.dt);
  Parameter::time++;

  param.DebugDumpAllParam();

  //main loop
  const PS::U32 atime = Parameter::time + Parameter::all_time - 1;
  for (; Parameter::time < atime; Parameter::time++) {
    drift_and_predict(system, param.dt, param.box_leng, param.ibox_leng);
    CHECK_VAL(system, DRIFT_KICK);

    dinfo.decomposeDomainAll(system);
    CHECK_VAL(system, DECOMPOSE);

    system.exchangeParticle(dinfo);
    CHECK_VAL(system, EXCHANGE);

    dens_tree.calcForceAllAndWriteBack(CalcDensity(), system, dinfo);
    force_tree.calcForceAllAndWriteBack(CalcForceEpEpDPD(), system, dinfo);
    CHECK_VAL(system, F_CALC_NB);

    kick(system, param.dt);
    CHECK_VAL(system, KICK);

    if (Parameter::time % Parameter::step_mac == 0) {
      observer.KineticTempera(system);
      observer.Pressure(system, bonded_vir, param.ibox_leng);
      observer.Diffusion(system, param.amp_num, param.sol_num);
      //observer.ConfigTempera();
    }

    if (Parameter::time % Parameter::step_mic == 0) {
      observer.Trajectory(system);
    }

    if (Parameter::time % Observer<PS::ParticleSystem<FPDPD> >::flush_freq == 0)
      observer.FlushAll();
  }//end of main loop

  //print configuration for restart
  observer.FinConfig(system);

  observer.CleanUp();

  const PS::S32 beg_time = std::atoi(argv[2]);
  calc_mean_pressure(cdir, beg_time, param.box_leng, param.box_leng.x * param.box_leng.z * 2.0 / param.init_amp_num);
  calc_mean_kintemp(cdir, beg_time, param.dt);

  PS::Finalize();
}
