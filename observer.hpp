#pragma once

#include <cstdio>
#include <vector>

template<class Tpsys>
class Observer {
  std::vector<FILE*> ptr_f;
  std::string cdir;

  void CleanUp() {
    for(auto it = ptr_f.begin(); it != ptr_f.end(); ++it)
      fclose(*it);
  }

  std::string type2fname(const PS::S32 type) const {
    std::string fname;
    switch(type){
    case KIN_TEMP:
      fname = "kin_temp";
      break;
    case CONF_TEMP:
      fname = "conf_temp";
      break;
    case PRESSURE:
      fname = "pressure";
      break;
    case DIFFUSION:
      fname = "diffus";
      break;
    case NUM_AMP:
      fname = "num_amp";
      break;
    case PART_CONFIG:
      fname = "traject";
    default:
      std::cerr << "Unknown type\n";
      std::cerr << "Error occurs at " __FILE__ << " " << __LINE__ << "." << std::endl;
      PS::Abort();
      break;
    }
    fname += ".txt";
    
    const std::string ret = cdir + "/" + fname;
    return ret;
  }

  //observer type
  enum {
    KIN_TEMP = 0,
    CONF_TEMP,
    PRESSURE,
    DIFFUSION,
    NUM_AMP,
    PART_CONFIG,
    
    NUM_FILES,
  };
  
public:
  Observer(const std::string& cdir_) {
    cdir = cdir_;
    ptr_f.resize(NUM_FILES, nullptr);
  }
  ~Observer() {
    CleanUp();
  }
  
  void Initialize() {
    for(PS::S32 i = 0; i < NUM_FILES; i++)
      ptr_f[i] = fopen(type2fname(i).c_str(), "w");
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
    fprintf(ptr_f[KIN_TEMP], "%.15g %.15g %.15g\n", kin_sum.x, kin_sum.y, kin_sum.z);
  }

  void ConfigTempera() {
    //not implemented yet
  }
  
  void Pressure(Tpsys& sys, PS::F64vec& ibox_leng) {
    //TODO: implement line/surface tension calculation part.
    const PS::S32 num_part = sys.getNumberOfParticleLocal();
    PS::F64vec press_sum(0.0, 0.0, 0.0);
    for(PS::S32 i = 0; i < num_part; i++) {
      press_sum += sys[i].press;
      press_sum.x += sys[i].vel.x * sys[i].vel.x;
      press_sum.y += sys[i].vel.y * sys[i].vel.y;
      press_sum.z += sys[i].vel.z * sys[i].vel.z;
    }
    press_sum *= ibox_leng.x * ibox_leng.y * ibox_leng.z;
    fprintf(ptr_f[PRESSURE], "%.15g %.15g %.15g\n", press_sum.x, press_sum.y, press_sum.z);
  }
  
  void Diffusion() {
    //not implemented yet
  }

  void Trajectory(const Tpsys& sys) {
    const PS::S32 num_part = sys.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < num_part; i++)
      sys[i].writeAscii(ptr_f[PART_CONFIG]);
  }
};
