#pragma once

#include <string>
#include <memory>
#include <map>
#include "ptcl_connector.hpp"
#include "particle_simulator.hpp"

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
  
  double est_grid_leng_ = 1.2, cutof_leng_    = 1.2;
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
    std::map<std::string, std::vector<std::string> > tag_val;
    Parameter::ReadTagValues(fin_run, tag_val);
    PS::F64vec box_leng;
    Parameter::Matching(&(box_leng[0]), std::string("box_leng"), tag_val, 3);
    std::cerr << "Read parameters.\n";

    ptr_connector.reset(new PtclConnector<PtclOrg, Ptcl>(box_leng.x, box_leng.y, box_leng.z));
    ptr_connector->Initialize(est_grid_leng_, cutof_leng_);
  }
  
  void SetSearchRadius(const double gl, const double cl) {
    assert(gl >= cl && "est_grid_leng_ should be greater than cutof_leng.");
    est_grid_leng_ = gl;
    cutof_leng_    = cl;
  }
  
  virtual void DoAnalysis() = 0;
};


#undef FILE_OPEN_WITH_CHECK
