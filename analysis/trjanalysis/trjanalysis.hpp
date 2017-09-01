#pragma once

#include <string>
#include <memory>
#include <map>
#include <functional>
#include "ptcl_connector.hpp"
#include "particle_simulator.hpp"
#include "user_defs.h"
#include "parameter.hpp"

#define FILE_OPEN_WITH_CHECK(f, fname)                        \
  do {                                                        \
    f.open(fname.c_str());                                    \
    if (f.fail()) {                                           \
      std::cerr << "Cannot open file: " << fname << ".\n";    \
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;	\
      std::exit(1);                                           \
    }                                                         \
  } while(false)

template<class PtclOrg, class Ptcl>
class TrjAnalysis {
protected:
  std::ifstream fin_trj, fin_run;
  const std::string cur_dir_;
  std::unique_ptr<PtclConnector<PtclOrg, Ptcl> > ptr_connector;

  // simulation info
  std::unique_ptr<Parameter> ptr_parameter;

  // ptcl data
  std::vector<Ptcl> ptcls_;
  std::vector<PtclOrg> ptcls_org_;

  // for neighbor search
  double est_grid_leng_ = 1.0, cutof_leng_ = 1.0;

  virtual void DebugDump()  {
    std::ofstream fout_trj("trj.xyz");

    fout_trj << ptcls_.size() << std::endl;
    fout_trj << "revised image\n";
    for (const auto& ptcl : ptcls_)
      fout_trj << "O " << ptcl.r << std::endl;
  }

  int GetLargestPatchId(const std::vector<int>& patch_id,
			const int num_patch) {
    std::vector<int> elem_in_each;
    for (int i = 0; i < num_patch; i++) {
      elem_in_each.push_back(std::count(patch_id.cbegin(), patch_id.cend(), i));
    }
    return std::distance(elem_in_each.begin(),
			 std::max_element(elem_in_each.begin(), elem_in_each.end()));
  }

  void GetNewBox(PS::F64vec& box_org, PS::F64vec& box_top) {
    PS::F64 min_x = std::numeric_limits<PS::F64>::max(), max_x = std::numeric_limits<PS::F64>::min();
    PS::F64 min_y = std::numeric_limits<PS::F64>::max(), max_y = std::numeric_limits<PS::F64>::min();
    PS::F64 min_z = std::numeric_limits<PS::F64>::max(), max_z = std::numeric_limits<PS::F64>::min();

    for (const auto& ptcl : ptcls_) {
      if (ptcl.r.x < min_x) min_x = ptcl.r.x;
      if (ptcl.r.x > max_x) max_x = ptcl.r.x;
      if (ptcl.r.y < min_y) min_y = ptcl.r.y;
      if (ptcl.r.y > max_y) max_y = ptcl.r.y;
      if (ptcl.r.z < min_z) min_z = ptcl.r.z;
      if (ptcl.r.z > max_z) max_z = ptcl.r.z;
    }

    box_org.x = min_x; box_org.y = min_y; box_org.z = min_z;
    box_top.x = max_x; box_top.y = max_y; box_top.z = max_z;

    const PS::F64vec offset(0.0);
    box_org -= offset;
    box_top += offset;
  }

  void ChangeCoordOrigin(const PS::F64vec box_org) {
    for (auto& ptcl : ptcls_) ptcl.r -= box_org;
  }

public:
  TrjAnalysis(const std::string cur_dir,
	      const char* trj_fname,
	      const char* param_fname) : cur_dir_(cur_dir) {
    const std::string fname_trj = cur_dir + trj_fname;
    const std::string fname_run = cur_dir + param_fname;
    FILE_OPEN_WITH_CHECK(fin_trj, fname_trj);
    FILE_OPEN_WITH_CHECK(fin_run, fname_run);
  }

  virtual ~TrjAnalysis() {}

  void Initialize() {
    ptr_parameter.reset(new Parameter(cur_dir_));
    ptr_parameter->Initialize();
    ptr_parameter->LoadParam();
    ptr_connector.reset(new PtclConnector<PtclOrg, Ptcl>);
  }

  void SetSearchRadius(const double gl, const double cl) {
    assert(gl >= cl && "est_grid_leng_ should be greater than cutof_leng.");
    est_grid_leng_ = gl;
    cutof_leng_    = cl;
  }

  bool ReadOneFrame(PS::U32& time,
		    const std::function<bool(const Ptcl&)>& is_target) {
    std::string tag, line;
    char buf[4];
    std::getline(fin_trj, line);
    if (fin_trj.eof()) return true;
    std::size_t num_ptcl = std::stoi(line);
    std::getline(fin_trj, line);
    sscanf(line.c_str(), "%s %u\n", buf, &time);

    ptcls_org_.resize(num_ptcl);
    for (std::size_t i = 0; i < num_ptcl; i++) {
      std::getline(fin_trj, line);
      ptcls_org_[i].readFromString(line.c_str());
    }

    // copy to buffer
    ptcls_.clear();
    Ptcl p;
    for (std::size_t i = 0; i < ptcls_org_.size(); i++) {
      p.CopyFromPtclOrg(ptcls_org_[i]);
      if (is_target(p)) ptcls_.push_back(p);
    }
    return false;
  }

  virtual void DoAnalysis() = 0;
};

#undef FILE_OPEN_WITH_CHECK
