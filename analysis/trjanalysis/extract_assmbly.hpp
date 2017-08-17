#pragma once
#include "trjanalysis.hpp"
#include "io_util.hpp"

class TrjAnalysisExtractAssmbly : public TrjAnalysis<FPDPD, Particle> {
  void WriteOneFrame(FILE* fp, const int num_patch) {
    const auto& ptch_id = ptr_connector->patch_id();
    const int tar_patch_id = GetLargestPatchId(ptch_id, num_patch);

    int num_tar_ptcl = 0;
    for (size_t i = 0; i < ptcls_.size(); i++) {
      if (ptch_id[i] == tar_patch_id) {
        num_tar_ptcl++;
      }
    }

    fprintf(fp, "%d\n", num_tar_ptcl);
    fprintf(fp, "molecular self-assembled structure is extracted.\n");
    for (size_t i = 0; i < ptcls_.size(); i++) {
      if (ptch_id[i] == tar_patch_id) {
        ptcls_org_[i].writeAscii(fp);
      }
    }
  }

public:
  TrjAnalysisExtractAssmbly(const std::string cur_dir,
                            const char* trj_fname,
                            const char* param_fname) : TrjAnalysis(cur_dir, trj_fname, param_fname) {}
  ~TrjAnalysisExtractAssmbly() override {}

  void DoAnalysis() override {
    SetSearchRadius(1.5, 1.5);
    ptr_connector->Initialize(est_grid_leng_, cutof_leng_, Parameter::box_leng);
    const std::string fname = cur_dir_ + "/extracted_mol_assembly.xyz";
    FILE* fp = io_util::xfopen(fname.c_str(), "w");

    PS::U32 time = 0;
    while (true) {
      if (ReadOneFrame(time, [] (const Particle& p) -> bool {return ((p.prop == Parameter::Hyphil) || (p.prop == Parameter::Hyphob));})) break;
      std::cout << "time = " << time << std::endl;
      ptr_connector->ResizeIfNeeded(ptcls_.size());

      const auto num_patch = gen_connected_image(ptr_connector.get(), ptcls_, Parameter::box_leng, false);
      std::cout << "# of patch is " << num_patch << std::endl;
      WriteOneFrame(fp, num_patch);
    }
    fclose(fp);
  }
};
