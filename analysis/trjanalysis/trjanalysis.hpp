#pragma once

#include <string>
#include <memory>
#include <map>
#include <functional>
#include "ptcl_connector.hpp"
#include "particle_simulator.hpp"
#include "user_defs.h"
#include "parameter.hpp"

#define FILE_OPEN_WITH_CHECK(f, fname)				\
  do {								\
    f.open(fname.c_str());					\
    if (f.fail()) {						\
      std::cerr << "Cannot open file: " << fname << ".\n";	\
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;	\
      std::exit(1);						\
    }								\
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
