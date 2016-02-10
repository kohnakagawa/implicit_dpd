#pragma once

#include <cstdio>
#include <vector>

#include <sstream>

template<class Tpsys>
class Observer {
  //observer type
  enum {
    KIN_TEMP = 0,
    CONF_TEMP,
    PRESSURE,
    DIFFUSION,
    NUM_AMP,
    PART_CONFIG,
    FIN_CONFIG,
    
    NUM_FILES,
  };

  FILE* ptr_f[NUM_FILES] = {nullptr};
  std::string cdir;

  void type2fname(const PS::S32 type, std::string& fname) {
    switch(type){
    case KIN_TEMP:
      fname = "kin_temp.txt";
      break;
    case CONF_TEMP:
      fname = "conf_temp.txt";
      break;
    case PRESSURE:
      fname = "pressure.txt";
      break;
    case DIFFUSION:
      fname = "diffus.txt";
      break;
    case NUM_AMP:
      fname = "num_amp.txt";
      break;
    case PART_CONFIG:
      fname = "traject.xyz";
      break;
    case FIN_CONFIG:
      fname = "fin_config.xyz";
      break;
    default:
      std::cerr << "Unknown type\n";
      std::cerr << "Error occurs at " __FILE__ << " " << __LINE__ << "." << std::endl;
      PS::Abort();
      break;
    }
    fname = cdir + "/" + fname;
  }
public:
  static constexpr PS::U32 flush_freq = 200;

  explicit Observer(const std::string cdir_) {
    cdir = cdir_;
  }
  ~Observer() {}
  
  void Initialize() {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    if(PS::Comm::getRank() == 0) {
#endif
      std::string fname;
      for(PS::U32 i = 0; i < NUM_FILES; i++) {
	type2fname(i, fname);
	ptr_f[i] = io_util::xfopen(fname.c_str(), "w");
      }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    }
#endif
  }
  
  void KineticTempera(const Tpsys& sys) {
    const PS::S32 num_loc = sys.getNumberOfParticleLocal();
    PS::F64vec kin_sum_loc(0.0, 0.0, 0.0);
    for(PS::S32 i = 0; i < num_loc; i++) {
      kin_sum_loc.x += sys[i].vel.x * sys[i].vel.x;
      kin_sum_loc.y += sys[i].vel.y * sys[i].vel.y;
      kin_sum_loc.z += sys[i].vel.z * sys[i].vel.z;
    }
    kin_sum_loc /= num_loc;

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    PS::F64vec kin_sum = PS::Comm::getSum(kin_sum_loc);
    kin_sum /= PS::Comm::getNumberOfProc();
    const PS::F64 Tmean = (kin_sum.x + kin_sum.y + kin_sum.z) / 3.0;
    if(PS::Comm::getRank() == 0)
      fprintf(ptr_f[KIN_TEMP], "%.15g %.15g %.15g %.15g\n", kin_sum.x, kin_sum.y, kin_sum.z, Tmean);

#else
    const PS::F64 Tmean_loc = (kin_sum_loc.x + kin_sum_loc.y + kin_sum_loc.z) / 3.0;
    fprintf(ptr_f[KIN_TEMP], "%.15g %.15g %.15g %.15g\n", kin_sum_loc.x, kin_sum_loc.y, kin_sum_loc.z, Tmean_loc);
#endif
  }

  void ConfigTempera() {
    //not implemented yet
  }
  
  void Pressure(const Tpsys& sys, const PS::F64vec& bonded_vir, const PS::F64vec& ibox_leng) {
    //NOTE: We do not use the action-reaction law when calculating the non-bonded interactions.
    //      Therefore, the virial should be multiplied by 0.5.
    const PS::S32 num_loc = sys.getNumberOfParticleLocal();
    PS::F64vec press_sum_loc(0.0, 0.0, 0.0);
    for(PS::S32 i = 0; i < num_loc; i++) {
      press_sum_loc += sys[i].press * 0.5;
      press_sum_loc.x += sys[i].vel.x * sys[i].vel.x;
      press_sum_loc.y += sys[i].vel.y * sys[i].vel.y;
      press_sum_loc.z += sys[i].vel.z * sys[i].vel.z;
    }
    press_sum_loc += bonded_vir * 0.5;
    
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    PS::F64vec press_sum = PS::Comm::getSum(press_sum_loc);
    press_sum *= ibox_leng.x * ibox_leng.y * ibox_leng.z;
    if(PS::Comm::getRank() == 0)
      fprintf(ptr_f[PRESSURE], "%.15g %.15g %.15g\n", press_sum.x, press_sum.y, press_sum.z);
#else
    press_sum_loc *= ibox_leng.x * ibox_leng.y * ibox_leng.z;
    fprintf(ptr_f[PRESSURE], "%.15g %.15g %.15g\n", press_sum_loc.x, press_sum_loc.y, press_sum_loc.z);
#endif
  }
  
  void Diffusion(const Tpsys& sys, const PS::U32 amp_num) {
    PS::F64 difsum_loc = 0.0;
    const PS::S32 num_loc = sys.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < num_loc; i++)
      difsum_loc += sys[i].delta_sumr * sys[i].delta_sumr;
    difsum_loc /= num_loc;

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
    //ASSUME: Molecular topology is constant during simulation.
    //        This routine is valid for constant amphiphile number simulation.
    assert(num_loc % Parameter::all_unit == 0);
    
    PS::F64 difsum_mol = 0.0;
    for(PS::U32 i = 0; i < amp_num; i++) {
      PS::F64vec dsumr_mol = 0.0;
      for(PS::U32 k = 0; k < Parameter::all_unit; k++)
	dsumr_mol += sys[Parameter::all_unit * i + k].delta_sumr;
      dsumr_mol /= Parameter::all_unit;
      difsum_mol += dsumr_mol * dsumr_mol;
    }
    difsum_mol /= amp_num;
    fprintf(ptr_f[DIFFUSION], "%.15g %.15g\n", difsum_loc, difsum_mol);
#elif defined PARTICLE_SIMULATOR_MPI_PARALLEL
    PS::F64 difsum = PS::Comm::getSum(difsum_loc);
    difsum /= PS::Comm::getNumberOfProc();
    if(PS::Comm::getRank() == 0)
      fprintf(ptr_f[DIFFUSION], "%.15g %.15g\n", difsum_loc, difsum);
#endif
  }

  void Trajectory(const Tpsys& sys) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    io_util::WriteXYZFormMPI(sys,
			     Parameter::time,
			     3 * sys.getNumberOfParticleGlobal(),
			     ptr_f[PART_CONFIG]); //TODO: getNumberOfParticleGlobal should be replaced.
#else
    io_util::WriteXYZForm(sys,
			  sys.getNumberOfParticleLocal(),
			  Parameter::time,
			  ptr_f[PART_CONFIG]);
#endif
  }

  void FinConfig(const Tpsys& sys) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    io_util::WriteXYZFormMPI(sys,
			     Parameter::time,
			     3 * sys.getNumberOfParticleGlobal(),
			     ptr_f[FIN_CONFIG]); //TODO: getNumberOfParticleGlobal should be replaced.
#else
    io_util::WriteXYZForm(sys,
			  sys.getNumberOfParticleLocal(),
			  Parameter::time,
			  ptr_f[FIN_CONFIG]);
#endif
  }

  void FlushAll() {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    if(PS::Comm::getRank() == 0) {
#endif
      for(PS::S32 i = 0; i < NUM_FILES; i++)
	fflush(ptr_f[i]);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    }
#endif 
  }

  void CleanUp() {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    if(PS::Comm::getRank() == 0) {
#endif
      for(PS::U32 i = 0; i < NUM_FILES; i++)
	fclose(ptr_f[i]);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    }
#endif 
    
  }
};
