#pragma once
#include "trjanalysis.hpp"
#include "io_util.hpp"

class TrjAnalysisExtractAssmbly : public TrjAnalysis<FPDPD, Particle> {
  void WriteOneFrame(FILE* fp, const int num_patch, const PS::U32 time) {
    const auto& ptch_id = ptr_connector->patch_id();
    const int tar_patch_id = GetLargestPatchId(ptch_id, num_patch);

    std::vector<FPDPD> ptcls_buffer;
    for (size_t i = 0; i < ptcls_org_.size(); i++) {
      if (ptch_id[i] == tar_patch_id) {
        ptcls_buffer.push_back(ptcls_org_[i]);
      }
    }
    std::sort(ptcls_buffer.begin(), ptcls_buffer.end(),
              [](const FPDPD& fp1, const FPDPD& fp2) -> bool {
                return (fp1.amp_id * Parameter::all_unit + fp1.unit) < (fp2.amp_id * Parameter::all_unit + fp2.unit);
              }
              );

    for (size_t i = 0; i < ptcls_buffer.size(); i++) {
      ptcls_buffer[i].id     = i;
      const size_t new_amp_id   = i / Parameter::all_unit;
      const size_t new_amp_unit = i % Parameter::all_unit;
      ptcls_buffer[i].amp_id = new_amp_id;
      ptcls_buffer[i].unit   = new_amp_unit;

      if (new_amp_unit < Parameter::head_unit) {
        assert(ptcls_buffer[i].prop == Parameter::Hyphil);
      } else {
        assert(ptcls_buffer[i].prop == Parameter::Hyphob);
      }
    }

    // NOTE: We do not consider solvent particles
    fprintf(fp, "%u\n", ptcls_buffer.size());
    fprintf(fp, "time %u\n", time);
    for (size_t i = 0; i < ptcls_buffer.size(); i++) {
      ptcls_buffer[i].writeAscii(fp);
    }
  }

public:
  TrjAnalysisExtractAssmbly(const std::string cur_dir,
                            const char* trj_fname,
                            const char* param_fname) : TrjAnalysis(cur_dir, trj_fname, param_fname) {}
  ~TrjAnalysisExtractAssmbly() override {}

  // ASSUME: there are no solvent particles.
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
      WriteOneFrame(fp, num_patch, time);
    }
    fclose(fp);
  }
};
