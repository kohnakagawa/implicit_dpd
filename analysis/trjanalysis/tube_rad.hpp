#pragma once
#include "trjanalysis.hpp"
#include "adjust_axis.hpp"

class TrjAnalysisTubeRad : public TrjAnalysis<FPDPD, Particle> {
  std::ofstream fout_;
  
  static constexpr int num_slice = 20;
  double slice_len_ = -1.0, tube_len_ = -1.0, min_z_ = -1.0, max_z_ = -1.0;
  std::array<PS::F64vec, num_slice> cm_in_slab_;
  std::array<double, num_slice> rad_in_slab_;
  std::array<int, num_slice> num_in_slab_;

  void CalcTubeLength(const int tar_patch_id,
		      const std::vector<int>& ptch_id) {
    min_z_ = std::numeric_limits<double>::max();
    max_z_ = std::numeric_limits<double>::min();
    for (size_t i = 0; i < ptcls_.size(); i++) {
      if (ptch_id[i] == tar_patch_id) {
	if (ptcls_[i].r.z < min_z_) min_z_ = ptcls_[i].r.z;
	if (ptcls_[i].r.z > max_z_) max_z_ = ptcls_[i].r.z;
      }
    }
    tube_len_ = max_z_ - min_z_;
    slice_len_ = tube_len_ / num_slice;
  }
  
  void CalcTubeRadius(const int tar_patch_id,
		      const std::vector<int>& ptch_id) {
    cm_in_slab_.fill(PS::F64vec(0.0, 0.0, 0.0));
    rad_in_slab_.fill(0.0); num_in_slab_.fill(0);
    
    for (size_t i = 0; i < ptcls_.size(); i++) {
      if (ptch_id[i] == tar_patch_id) {
	int zid = static_cast<int>((ptcls_[i].r.z - min_z_) / slice_len_);
	if (zid == num_slice) zid--;
	cm_in_slab_[zid] += ptcls_[i].r;
	num_in_slab_[zid]++;
      }
    }

    for (int i = 0; i < num_slice; i++) {
      cm_in_slab_[i] /= num_in_slab_[i];
    }
    
    for (size_t i = 0; i < ptcls_.size(); i++) {
      if (ptch_id[i] == tar_patch_id) {
	int zid = static_cast<int>((ptcls_[i].r.z - min_z_) / slice_len_);
	if (zid == num_slice) zid--;
	const auto cm2ptcl = ptcls_[i].r - cm_in_slab_[zid];
	rad_in_slab_[zid] += std::sqrt(cm2ptcl.x * cm2ptcl.x + cm2ptcl.y * cm2ptcl.y);
      }
    }
    
    for (int i = 0; i < num_slice; i++) {
      rad_in_slab_[i] /= num_in_slab_[i];
    }
  }

  void WriteTubeRadDistrib() {
    for (int i = 0; i < num_slice; i++) {
      fout_ << (i + 0.5) * slice_len_ << " " << rad_in_slab_[i] << std::endl;
    }
    fout_ << std::endl;
  }
  
public:
  TrjAnalysisTubeRad(const std::string cur_dir,
		     const char* trj_fname,
		     const char* param_fname) : TrjAnalysis(cur_dir, trj_fname, param_fname) {}
  
  ~TrjAnalysisTubeRad() override {}
  
  void DoAnalysis() override {
    SetSearchRadius(1.2, 1.2);
    ptr_connector->Initialize(est_grid_leng_, cutof_leng_, Parameter::box_leng);
    
    const std::string fname = cur_dir_ + "/slice_rad.txt";
    fout_.open(fname.c_str());
   
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
      
      CalcTubeLength(tar_patch_id, ptch_id);
      CalcTubeRadius(tar_patch_id, ptch_id);
      WriteTubeRadDistrib();

      ptr_connector->ClearForNextStep();
    }

    DebugDump();
  }
};
