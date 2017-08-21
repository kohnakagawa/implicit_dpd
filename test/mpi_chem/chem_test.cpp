#include <iostream>
#include <iomanip>
#include <sstream>
#include "particle_simulator.hpp"
#include "io_util.hpp"
#include "parameter.hpp"
#include "f_calculator.hpp"
#include "chemmanager.hpp"
#include "observer.hpp"

constexpr char Parameter::atom_type[21];

PS::F64vec Parameter::box_leng, Parameter::ibox_leng;
PS::U32 Parameter::time;
PS::U32 Parameter::all_time, Parameter::step_mic, Parameter::step_mac;

template <class Tpsys>
void clear_acc(Tpsys& system) {
  const PS::S32 numptcls = system.getNumberOfParticleLocal();
  for (PS::S32 i = 0; i < numptcls; i++) {
    system[i].acc = 0.0;
  }
}

template <class Tpsys>
void print_added_ptcls(Tpsys& system,
                       const PS::S32 init_ptcl_num,
                       const std::string cdir) {
  const PS::S32 numptcls = system.getNumberOfParticleLocal();

  const PS::S32 rank = PS::Comm::getRank();
  std::stringstream ss;
  ss << cdir << "/rank" << rank << ".xyz";
  std::ofstream fout(ss.str().c_str());

  fout << std::setprecision(15);
  for (PS::S32 i = 0; i < numptcls; i++) {
    if (system[i].id >= init_ptcl_num) {
      fout << "O " << system[i].pos << std::endl;
    }
  }
}

int main(int argc, char* argv[]) {
  PS::Initialize(argc, argv);

  if (argc != 2) {
    std::cerr << "argv[1] = dirname" << std::endl;
    PS::Abort();
  }

  const std::string cdir = argv[1];
  Parameter param(cdir);
  param.Initialize();
  if (PS::Comm::getRank() == 0) param.LoadParam();

  PS::ParticleSystem<FPDPD> system;
  system.initialize();
  if (PS::Comm::getRank() == 0) {
    system.setNumberOfParticleLocal(param.init_amp_num * Parameter::all_unit + param.sol_num);
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

  PS::TreeForForceShort<RESULT::Density, EPI::Density, EPJ::Density>::Gather dens_tree;
  PS::TreeForForceShort<RESULT::ForceDPD, EPI::DPD, EPJ::DPD>::Gather force_tree;
  dens_tree.initialize(3 * system.getNumberOfParticleGlobal() );
  force_tree.initialize(3 * system.getNumberOfParticleGlobal() );
  dens_tree.calcForceAllAndWriteBack(CalcDensity(), system, dinfo);
  force_tree.calcForceAllAndWriteBack(CalcForceEpEpDPD(), system, dinfo);

  clear_acc(system);

  PS::F64vec bonded_vir(0.0, 0.0, 0.0);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  const PS::U32 est_loc_amp = param.init_amp_num / PS::Comm::getNumberOfProc() + 100;
  ForceBondedMPI<PS::ParticleSystem<FPDPD>, EPJ::DPD> fbonded(est_loc_amp);
  fbonded.CalcListedForce(system, force_tree.epj_org(), bonded_vir, param.core_amp_id());
#else
  ForceBonded<PS::ParticleSystem<FPDPD> > fbonded(system, Parameter::all_unit * param.init_amp_num, param.init_amp_num, param.core_amp_id());
  fbonded.CalcListedForce(system, bonded_vir);
#endif

  // main test part
  const PS::S32 init_ptcl_num = system.getNumberOfParticleGlobal();
  FILE* fp = nullptr;
  if (PS::Comm::getRank() == 0) {
    const std::string fname = cdir + "/increased_num.dat";
    fp = fopen(fname.c_str(), "w");
  }

  const PS::S32 est_max_ptcl_num = param.sol_num + static_cast<PS::S32>(param.max_amp_num * Parameter::all_unit * 1.3);
  ChemManager<FPDPD> chemmanag(1234);
  system.ExpandParticleBuffer(est_max_ptcl_num);
#ifndef PARTICLE_SIMULATOR_MPI_PARALLEL
  fbonded.ExpandTopolBuffer(param.max_amp_num * Parameter::all_unit);
#endif

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  (void) chemmanag.RandomChemEvent(system, force_tree.epj_org(),
                                   fbonded.loc_topol_cmpl(), fbonded.loc_topol_imcmpl(),
                                   fbonded.cmplt_amp(), fbonded.imcmplt_amp(),
                                   param);
#else
  (void) chemmanag.RandomChemEvent(system, fbonded.glob_topol(), param);
#endif

  const PS::S32 increased_num = system.getNumberOfParticleGlobal() - init_ptcl_num;
  if (PS::Comm::getRank() == 0) fprintf(fp, "%d\n", increased_num);

  print_added_ptcls(system, init_ptcl_num, cdir);

  Observer<PS::ParticleSystem<FPDPD> > observer(cdir);
  observer.Initialize();
  observer.Trajectory(system);

  PS::Finalize();
}
