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
    std::string fname;
    for(PS::U32 i = 0; i < NUM_FILES; i++) {
      type2fname(i, fname);
      ptr_f[i] = io_util::xfopen(fname.c_str(), "w");
    }
  }
  
  void KineticTempera(Tpsys& sys) {
    const PS::S32 num_part = sys.getNumberOfParticleLocal();
    PS::F64vec kin_sum(0.0, 0.0, 0.0);
    for(PS::S32 i = 0; i < num_part; i++) {
      kin_sum.x += sys[i].vel.x * sys[i].vel.x;
      kin_sum.y += sys[i].vel.y * sys[i].vel.y;
      kin_sum.z += sys[i].vel.z * sys[i].vel.z;
    }
    kin_sum /= num_part;
    const PS::F64 Tmean = (kin_sum.x + kin_sum.y + kin_sum.z) / 3.0;
    fprintf(ptr_f[KIN_TEMP], "%.15g %.15g %.15g %.15g\n", kin_sum.x, kin_sum.y, kin_sum.z, Tmean);
  }

  void ConfigTempera() {
    //not implemented yet
  }
  
  void Pressure(Tpsys& sys, PS::F64vec& bonded_vir, PS::F64vec& ibox_leng) {
    //NOTE: We do not use the action-reaction law when calculating the non-bonded interactions.
    //      Therefore, the virial should be multiplied by 0.5.
    //TODO: implement line/surface tension calculation part.
    const PS::S32 num_part = sys.getNumberOfParticleLocal();
    PS::F64vec press_sum(0.0, 0.0, 0.0);
    for(PS::S32 i = 0; i < num_part; i++) {
      press_sum += sys[i].press * 0.5;
      press_sum.x += sys[i].vel.x * sys[i].vel.x;
      press_sum.y += sys[i].vel.y * sys[i].vel.y;
      press_sum.z += sys[i].vel.z * sys[i].vel.z;
    }
    press_sum += bonded_vir * 0.5;
    
    press_sum *= ibox_leng.x * ibox_leng.y * ibox_leng.z;
    fprintf(ptr_f[PRESSURE], "%.15g %.15g %.15g\n", press_sum.x, press_sum.y, press_sum.z);
  }
  
  void Diffusion(const Tpsys& sys, const PS::S32 amp_num) {
    PS::F64 difsum = 0.0;
    const PS::S32 num = sys.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < num; i++)
      difsum += sys[i].delta_sumr * sys[i].delta_sumr;
    difsum /= num;

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
    //ASSUME: Molecular topology is constant during simulation.
    //        This routine is valid for constant amphiphile number simulation.
    assert(num % Parameter::all_unit == 0);
    
    PS::F64 difsum_mol = 0.0;
    for(PS::S32 i = 0; i < amp_num; i++) {
      PS::F64vec dsumr_mol = 0.0;
      for(PS::S32 k = 0; k < Parameter::all_unit; k++)
	dsumr_mol += sys[Parameter::all_unit * i + k].delta_sumr;
      dsumr_mol /= Parameter::all_unit;
      difsum_mol += dsumr_mol * dsumr_mol;
    }
    difsum_mol /= amp_num;
    fprintf(ptr_f[DIFFUSION], "%.15g %.15g\n", difsum, difsum_mol);
#else
    fprintf(ptr_f[DIFFUSION], "%.15g\n", difsum);
#endif
  }

  void Trajectory(const Tpsys& sys) {
    io_util::WriteXYZForm(sys,
			  sys.getNumberOfParticleLocal(),
			  Parameter::time,
			  ptr_f[PART_CONFIG]);
  }

  void FinConfig(const Tpsys& sys) {
    io_util::WriteXYZForm(sys,
			  sys.getNumberOfParticleLocal(),
			  Parameter::time,
			  ptr_f[FIN_CONFIG]);
  }

  void FlushAll() {
    for(PS::S32 i = 0; i < NUM_FILES; i++)
      fflush(ptr_f[i]);
  }

  void CleanUp() {
    for(PS::U32 i = 0; i < NUM_FILES; i++)
      fclose(ptr_f[i]);
  }
};
