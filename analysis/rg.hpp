#pragma once
#include "trjanalysis.hpp"

class TrjAnalysisRg : public TrjAnalysis<FPDPD, Particle> {
  double CalcRg(const std::vector<Particle>& ptcls) {
    PS::F64vec cm_pos = 0.0;
    for (const auto& ptcl : ptcls) cm_pos += ptcl.r;
    cm_pos /= ptcls.size();

    double rg = 0.0;
    for (const auto& ptcl : ptcls) {
      const auto cm_to_pos = ptcl.r - cm_pos;
      rg += cm_to_pos * cm_to_pos;
    }
    rg = std::sqrt(rg / ptcls.size());
    return rg;
  }
  
public:
  TrjAnalysisRg(const std::string cur_dir,
		const char* trj_fname,
		const char* param_fname) : TrjAnalysis(cur_dir, trj_fname, param_fname) {}
  ~TrjAnalysisRg() override {}
  
  void DoAnalysis() override {
    const std::string name_fout = cur_dir_ + "/rg.txt";
    std::ofstream fout(name_fout.c_str());
    
    size_t time = 0;
    while (true) {
      if (ptr_connector->ReadOneFrame(fin_trj, time)) break;
      std::cout << "time = " << time << std::endl;
      const auto num_patch = gen_connected_image(ptr_connector.get(), false);
      std::cout << "# of patch is " << num_patch << std::endl;
      fout << time << " " << CalcRg(ptr_connector->ptcls());
      ptr_connector->ClearForNextStep();
    }
  }
};
