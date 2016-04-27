#pragma once
#include "trjanalysis.hpp"

class TrjAnalysisRg : public TrjAnalysis<FPDPD, Particle> {
  double CalcRg(const int num_patch) {
    const auto& ptch_id = ptr_connector->patch_id();

    const int tar_patch_id = GetLargestPatchId(ptch_id, num_patch);

    int num_tar_ptcl = 0;
    PS::F64vec cm_pos = 0.0;
    for (size_t i = 0; i < ptcls_.size(); i++) {
      if (ptch_id[i] == tar_patch_id) {
	cm_pos += ptcls_[i].r;
	num_tar_ptcl++;
      }
    }
    cm_pos /= num_tar_ptcl;

    double rg = 0.0;
    for (size_t i = 0; i < ptcls_.size(); i++) {
      if (ptch_id[i] == tar_patch_id) {
	const auto cm_to_pos = ptcls_[i].r - cm_pos;
	rg += cm_to_pos * cm_to_pos;
      }
    }
    rg = std::sqrt(rg / num_tar_ptcl);
    return rg;
  }
  
public:
  TrjAnalysisRg(const std::string cur_dir,
		const char* trj_fname,
		const char* param_fname) : TrjAnalysis(cur_dir, trj_fname, param_fname) {}
  ~TrjAnalysisRg() override {}
  
  void DoAnalysis() override {
    SetSearchRadius(1.5, 1.5);
    ptr_connector->Initialize(est_grid_leng_, cutof_leng_, Parameter::box_leng);    
    
    const std::string name_fout = cur_dir_ + "/rg.txt";
    std::ofstream fout_rg(name_fout.c_str());

    PS::U32 time = 0;
    while (true) {
      if (ReadOneFrame(time, [] (const Particle& p) -> bool {return ((p.prop == Parameter::Hyphil) || (p.prop == Parameter::Hyphob));})) break;
      std::cout << "time = " << time << std::endl;
      ptr_connector->ResizeIfNeeded(ptcls_.size());

      const auto num_patch = gen_connected_image(ptr_connector.get(), ptcls_, Parameter::box_leng, false);
      std::cout << "# of patch is " << num_patch << std::endl;
      
      fout_rg << time << " " << CalcRg(num_patch) << std::endl;
      
      // DebugDump();

      ptr_connector->ClearForNextStep();
    }
  }
};
