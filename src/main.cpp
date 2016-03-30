#include <iostream>
#include <sys/time.h>
#include "particle_simulator.hpp"
#include "user_defs.h"
#include "parameter.hpp"

#ifdef ENABLE_GPU_CUDA
#warning "will use gpu."
#include "ptcl_class.hpp"
#include "bdf_calculator.hpp"
#include "f_calculator_gpu.cuh"
#include "device_inform.cuh"
#else
#include "f_calculator.hpp"
#endif

#ifdef CHEM_MODE
#include "chemmanager.hpp"
#endif

#include "observer.hpp"
#include "driftkick.hpp"

#if defined (PARTICLE_SIMULATOR_MPI_PARALLEL) && defined (CALC_HEIGHT)
#error "MPI version of CALC_HEIGHT is not supported!"
#endif

#ifdef CALC_HEIGHT
#warning "will calculate the membrane height. This calculation causes performance loss."
#endif

#ifdef CHEM_MODE
#warning "Chemical reaction occurs."
#endif

constexpr char Parameter::atom_type[21];

PS::F64vec Parameter::box_leng, Parameter::ibox_leng;
PS::U32 Parameter::time;
PS::U32 Parameter::all_time, Parameter::step_mic, Parameter::step_mac;

static_assert(Parameter::head_unit == 1, "head_unit should be 1.");
static_assert(Parameter::tail_unit == 3, "tail_unit should be 3.");
static_assert(Parameter::bond_leng != 0.0, "bond_leng should be not 0.0.");
static_assert(Parameter::Reo < 3.0, "Reo should be less than 3.0.");

namespace {
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

  // helper functions for observer
  inline void do_observe_macro(Observer<PS::ParticleSystem<FPDPD> >& observer,
			       const PS::ParticleSystem<FPDPD>& system,
			       const Parameter& param,
			       const PS::F64vec& bonded_vir) {
    observer.KineticTempera(system);
    observer.Pressure(system, bonded_vir, param.ibox_leng);
    observer.Diffusion(system, param.amp_num, param.sol_num);
    
#ifdef CALC_HEIGHT
    observer.MembHeight(system, param.box_leng);
#endif

#ifdef CHEM_MODE
    observer.NumAmp(param.amp_num);
#endif
  }
  inline void do_observe_micro(Observer<PS::ParticleSystem<FPDPD> >& observer,
			       const PS::ParticleSystem<FPDPD>& system) {
    observer.Trajectory(system);
  }
} // end of anonymous namespace

int main(int argc, char *argv[]) {
  timer_start();

  PS::Initialize(argc, argv);
  
#ifdef ENABLE_GPU_CUDA
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  // Multi GPU case
  const PS::S32 dev_id = PS::Comm::getRank() % NUM_GPU_IN_ONE_NODE;
  checkCudaErrors(cudaSetDevice(dev_id));
#else
  // Single GPU case
  if (argc == 3) {
    const PS::S32 dev_id = std::atoi(argv[2]);
    checkCudaErrors(cudaSetDevice(dev_id));
    print_device_inform(dev_id);
  } else {
    std::cerr << "gpu device id is not specified.\n";
    PS::Abort();
  }
#endif
#endif

  // Initialize run input parameter
  const std::string cdir = argv[1];
  Parameter param(cdir);
  param.Initialize();
  if (PS::Comm::getRank() == 0) param.LoadParam();
  
  // Read restart configuration
  PS::ParticleSystem<FPDPD> system;
  system.initialize();
  if (PS::Comm::getRank() == 0) {
    system.setNumberOfParticleLocal(param.init_amp_num * Parameter::all_unit + param.sol_num);
    Parameter::time = param.LoadParticleConfig(system);
  } else {
    system.setNumberOfParticleLocal(0);
  }
  
#ifdef CHEM_MODE
  // Calc core id if needed.
  if (PS::Comm::getRank() == 0) param.CalcCorePtclId(system);
#endif

  // Share parameter data with other processes.
  param.ShareDataWithOtherProc();
  param.CheckLoaded();
  param.CheckParticleConfigIsValid(system);

  // Distribute particle data with other processes.
  PS::DomainInfo dinfo;
  const PS::F64 coef_ema = 0.3;
  dinfo.initialize(coef_ema);
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), Parameter::box_leng);
  dinfo.decomposeDomainAll(system);
  system.exchangeParticle(dinfo);

  // Initial step
  drift_and_predict(system, param.dt, param.box_leng, param.ibox_leng);
  dinfo.decomposeDomain();
  system.exchangeParticle(dinfo);
  PS::TreeForForceShort<RESULT::Density, EPI::Density, EPJ::Density>::Gather dens_tree;
  PS::TreeForForceShort<RESULT::ForceDPD, EPI::DPD, EPJ::DPD>::Gather force_tree;
  dens_tree.initialize(5 * system.getNumberOfParticleGlobal() );
  force_tree.initialize(5 * system.getNumberOfParticleGlobal() );
  
#ifdef ENABLE_GPU_CUDA
  const PS::S32 n_walk_limit = 800;
  const PS::S32 tag_max = 1;
  dens_tree.calcForceAllAndWriteBackMultiWalk(DispatchKernel<Policy::Density, EPI::Density, EPJ::Density>,
					      RetrieveKernel<Policy::Density, RESULT::Density>,
					      tag_max,
					      system,
					      dinfo,
					      n_walk_limit);
					      
  force_tree.calcForceAllAndWriteBackMultiWalk(DispatchKernel<Policy::Force, EPI::DPD, EPJ::DPD>,
					       RetrieveKernel<Policy::Force, RESULT::ForceDPD>,
					       tag_max,
					       system,
					       dinfo,
					       n_walk_limit);
#else
  // CalcForceEpEpDPD::m_seed = static_cast<PS::U32>(time(nullptr));
  CalcForceEpEpDPD::m_seed = 100;
  dens_tree.calcForceAllAndWriteBack(CalcDensity(), system, dinfo);
  force_tree.calcForceAllAndWriteBack(CalcForceEpEpDPD(), system, dinfo);
#endif

  PS::F64vec bonded_vir(0.0, 0.0, 0.0);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  const PS::U32 est_loc_amp = param.init_amp_num / PS::Comm::getNumberOfProc() + 100;
  ForceBondedMPI<PS::ParticleSystem<FPDPD>, EPJ::DPD> fbonded(est_loc_amp);
  fbonded.CalcListedForce(system, force_tree.epj_org(), bonded_vir);
#else
  ForceBonded<PS::ParticleSystem<FPDPD> > fbonded(system, Parameter::all_unit * param.init_amp_num);
  fbonded.CalcListedForce(system, bonded_vir);
#endif
  
  Observer<PS::ParticleSystem<FPDPD> > observer(cdir);
  observer.Initialize();
  do_observe_macro(observer, system, param, bonded_vir);
  do_observe_micro(observer, system);
#ifdef CHEM_MODE
#ifndef PARTICLE_SIMULATOR_MPI_PARALLEL
      observer.MembNormalVect(system, chemmanag.h2t_vecs(), chemmanag.core_poss_h(), param.core_ptcl_id);
#endif
#endif

#ifdef CHEM_MODE
  // const PS::U32 seed = 123;
  const PS::U32 seed = static_cast<PS::U32>(time(nullptr));
  ChemManager<FPDPD> chemmanag(seed);
  system.ExpandParticleBuffer(param.max_amp_num * Parameter::all_unit);
  
#ifndef PARTICLE_SIMULATOR_MPI_PARALLEL
  fbonded.ExpandTopolBuffer(param.max_amp_num * Parameter::all_unit);
#endif
  
#endif

  kick(system, param.dt);
  Parameter::time++;
  // End of initial step.
  
  // Main MD loop
  const PS::U32 atime = Parameter::time + Parameter::all_time - 1;
  for (; Parameter::time < atime; Parameter::time++) {
    drift_and_predict(system, param.dt, param.box_leng, param.ibox_leng);
    
    if (Parameter::time % Parameter::decom_freq == 0) dinfo.decomposeDomainAll(system);

    system.exchangeParticle(dinfo);

#ifdef ENABLE_GPU_CUDA
    dens_tree.calcForceAllAndWriteBackMultiWalk(DispatchKernel<Policy::Density, EPI::Density, EPJ::Density>,
					        RetrieveKernel<Policy::Density, RESULT::Density >,
					        tag_max,
					        system,
					        dinfo,
					        n_walk_limit);
					      
    force_tree.calcForceAllAndWriteBackMultiWalk(DispatchKernel<Policy::Force, EPI::DPD, EPJ::DPD>,
						 RetrieveKernel<Policy::Force, RESULT::ForceDPD>,
					         tag_max,
					         system,
					         dinfo,
					         n_walk_limit);
#else
    dens_tree.calcForceAllAndWriteBack(CalcDensity(), system, dinfo);
    force_tree.calcForceAllAndWriteBack(CalcForceEpEpDPD(), system, dinfo);
#endif

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    fbonded.CalcListedForce(system, force_tree.epj_org(), bonded_vir);
#else
    fbonded.CalcListedForce(system, bonded_vir);
#endif
    
    kick(system, param.dt);

#ifdef CHEM_MODE
    if (Parameter::time >= Parameter::beg_chem) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
      if (!chemmanag.RandomChemEvent(system, force_tree.epj_org(),
				     fbonded.loc_topol_cmpl(), fbonded.loc_topol_imcmpl(),
				     fbonded.cmplt_amp(), fbonded.imcmplt_amp(),
				     param)) break;
#else
      if (!chemmanag.RandomChemEvent(system, fbonded.glob_topol, param)) break;
#endif
    }
#endif

    if (Parameter::time % Parameter::step_mac == 0) {
      do_observe_macro(observer, system, param, bonded_vir);
#ifdef CHEM_MODE
#ifndef PARTICLE_SIMULATOR_MPI_PARALLEL
      observer.MembNormalVect(system, chemmanag.h2t_vecs(), chemmanag.core_poss_h(), param.core_ptcl_id);
#endif
#endif
    }

    if (Parameter::time % Parameter::step_mic == 0)
      do_observe_micro(observer, system);
  }
  // end of Main MD loop
  timer_stop();
  if (PS::Comm::getRank() == 0)
    show_duration();

  // print configuration for restart
  observer.FinConfig(system);

  // cleanup
  observer.CleanUp();
  param.DumpAllParam();
#ifdef ENABLE_GPU_CUDA
  clean_up_gpu();
#endif
  PS::Finalize();
}
