#pragma once

#include <cstdio>
#include <vector>
#include <string>
#include <map>

class Observer {
  std::map<std::string, FILE*> ptr_f;
  std::string cdir;

  void CleanUp() {
    for(auto it = ptr_f.begin(); it != ptr_f.end(); ++it)
      fclose(it->second);
  }
public:
  Observer(std::string cdir_) {
    cdir = cdir_;
  }
  ~Observer() {
    CleanUp();
  }
  
  void Initialize() {
    
    const std::string fname = cdir + "/kinetic.txt";
    ptr_f["kT"] = fopen(fname.c_str(), "w");
  }
  
  void KineticTempera() {
    
  }

  void ConfigTempera() {
    
  }
  
  void Pressure() {
    
  }
  
  void Diffusion() {
    
  }
};
