#pragma once
#include "trjanalysis.hpp"

class TrjAnalysisAde : public TrjAnalysis<FPDPD, Particle> {
  std::vector<int> CalcInOutElem(const std::vector<int>& patch_id,
				 const int num_patch) {
    std::vector<int> elem_in_each;
    for (int i = 0; i < num_patch; i++) {
      elem_in_each.push_back(std::count(patch_id.cbegin(), patch_id.cend(), i));
    }
    std::sort(elem_in_each.begin(), elem_in_each.end(), std::greater<int>());
    return elem_in_each;
  }

public:
  TrjAnalysisAde(const std::string cur_dir,
		 const char* trj_fname,
		 const char* param_fname) : TrjAnalysis(cur_dir, trj_fname, param_fname) {}
  ~TrjAnalysisAde() override {}

  void DoAnalysis() override {
    SetSearchRadius(1.1, 1.1);
    ptr_connector->Initialize(est_grid_leng_, cutof_leng_, Parameter::box_leng);
    
    const std::string name_fout = cur_dir_ + "/elem_in_each.txt";
    std::ofstream fout(name_fout.c_str());
    fout << "#time\tnum0\tnum1\thyphil_all\n";
    
    size_t time = 0;
    while (true) {
      if (ReadOneFrame(time, [] (const Particle& p) -> bool { return (p.prop == Parameter::Hyphil); })) break;
      std::cout << "time = " << time << std::endl;
      ptr_connector->ResizeIfNeeded(ptcls_.size());
      
      const auto num_patch = gen_connected_image(ptr_connector.get(), ptcls_, Parameter::box_leng, true);
      std::cout << "# of patch is " << num_patch << std::endl;
      
      const auto elem_in_each = CalcInOutElem(ptr_connector->patch_id(), num_patch);
      fout << time;
      for (const auto e : elem_in_each) fout << "\t" << e;
      fout << std::endl;
      
      ptr_connector->ClearForNextStep();
    }
  }
};
