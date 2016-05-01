#pragma once
#include "trjanalysis.hpp"
#include "adjust_axis.hpp"

class TrjAnalysisDensity : public TrjAnalysis<FPDPD, Particle> {
  std::ofstream fout_;
  static constexpr PS::F64 grid_len_est = 2.0;
  PS::F64vec grid_len_;
  std::array<int, 3> grid_n_;
  std::vector<PS::F64> num_in_box_;
  int nframes_ = 0;

  int GetHash(PS::F64vec r) {
    int id[] = {
      static_cast<int>(r.x / grid_len_.x),
      static_cast<int>(r.y / grid_len_.y),
      static_cast<int>(r.z / grid_len_.z),
    };
    for (auto i = 0; i < 3; i++) if (id[i] == grid_n_[i]) id[i]--;
    const auto ret = id[0] + grid_n_[0] * (id[1] + grid_n_[1] * id[2]);
    return ret;
  }

  void AddGridDensity() {
    for (const auto& ptcl : ptcls_) num_in_box_[GetHash(ptcl.r)]++;
    nframes_++;
  }

  void WriteVTKForm(const std::string fname) {
    std::ofstream fout(fname.c_str());
    fout << "# vtk DataFile Version 3.0.\n";
    fout << "vtk output\n";
    fout << "ASCII\n";
    fout << "DATASET STRUCTURED_POINTS\n";
    fout << "DIMENSIONS " << grid_n_[0] << " " << grid_n_[1] << " " << grid_n_[2] << std::endl;
    fout << "ORIGIN 0.0 0.0 0.0\n";
    fout << "SPACING " << grid_len_[0] << " " << grid_len_[1] << " " << grid_len_[2] << std::endl;
    fout << "POINT_DATA " << grid_n_[0] * grid_n_[1] * grid_n_[2] << std::endl;
    fout << "SCALARS scalars float\n";
    fout << "LOOKUP_TABLE default\n";
    for (const auto data : num_in_box_) fout << data << std::endl;
  }
  
public:
  TrjAnalysisDensity(const std::string cur_dir,
		     const char* trj_fname,
		     const char* param_fname) : TrjAnalysis(cur_dir, trj_fname, param_fname) {}
  
  ~TrjAnalysisDensity() override {}
  
  void DoAnalysis() override {
    SetSearchRadius(1.2, 1.2);
    ptr_connector->Initialize(est_grid_leng_, cutof_leng_, Parameter::box_leng);
    
    PS::U32 time = 0;
    AxisAdjuster<Particle> aadjuster;
    while (true) {
      if (ReadOneFrame(time, [] (const Particle& p) -> bool { return (p.prop == Parameter::Hyphil) || (p.prop == Parameter::Hyphob); })) break;
      std::cout << "time = " << time << std::endl;
      ptr_connector->ResizeIfNeeded(ptcls_.size());
      const auto num_patch = gen_connected_image(ptr_connector.get(), ptcls_, Parameter::box_leng, false);
      std::cout << "# of patch is " << num_patch << std::endl;

      const auto& ptch_id = ptr_connector->patch_id();
      const auto tar_patch_id = GetLargestPatchId(ptch_id, num_patch);
      aadjuster.CreateMomInertia(ptcls_, ptch_id, tar_patch_id);
      aadjuster.MakeNewBaseVector(ptcls_);
      aadjuster.ChangeBaseVector(ptcls_);
      
      PS::F64vec box_org(0.0, 0.0, 0.0), box_top(0.0, 0.0, 0.0);
      GetNewBox(box_org, box_top);
      ChangeCoordOrigin(box_org);
      const auto box_leng = box_top - box_org;
      for (auto i = 0; i < 3; i++) {
	grid_n_[i] = static_cast<int>(box_leng[i] / grid_len_est);
	grid_len_[i] = box_leng[i] / grid_n_[i];
      }
      const auto grid_n_all = grid_n_[0] * grid_n_[1] * grid_n_[2];
      num_in_box_.resize(grid_n_all, 0.0);
      AddGridDensity();
    }
    
    for (const auto g : grid_n_) std::cout << g << std::endl;
    
    const std::string fname = cur_dir_ + "/density_map.vtk";
    WriteVTKForm(fname);
  }
};
