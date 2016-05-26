#pragma once
#include "trjanalysis.hpp"

class TrjAnalysisReo : public TrjAnalysis<FPDPD, Particle> {
public:
  TrjAnalysisReo(const std::string cur_dir,
		 const char* trj_fname,
		 const char* param_fname) : TrjAnalysis(cur_dir, trj_fname, param_fname) {}
  
  ~TrjAnalysisReo() override {}
  
  void DoAnalysis() override {
    std::string fname = cur_dir_ + "/Reo.txt";
    std::ofstream fout;
    fout.open(fname.c_str());
   
    PS::U32 time = 0;
    while (true) {
      if (ReadOneFrame(time, [] (const Particle& p) -> bool { return (p.prop == Parameter::Hyphil) || (p.prop == Parameter::Hyphob); })) break;
      std::cout << "time = " << time << std::endl;
      
      double Reo2 = 0.0;
      const int amp_num = ptcls_.size() / Parameter::all_unit;
      for (int i = 0; i < amp_num; i++) {
	const auto r_h = ptcls_[i * Parameter::all_unit];
	const auto r_e = ptcls_[(i + 1) * Parameter::all_unit - 1];
	
	auto dr = r_e.r - r_h.r;
	dr[0] -= Parameter::box_leng[0] * std::round(dr[0] / Parameter::box_leng[0]);
	dr[1] -= Parameter::box_leng[1] * std::round(dr[1] / Parameter::box_leng[1]);
	dr[2] -= Parameter::box_leng[2] * std::round(dr[2] / Parameter::box_leng[2]);
	
	if (r_e.amp_id != r_h.amp_id) {
	  std::cerr << "Please check molecular topology.\n";
	  std::exit(1);
	}
	
	Reo2 += dr * dr;
	std::cout << std::sqrt(dr * dr) << std::endl;
	
      }
      Reo2 /= amp_num;
      
      fout << std::sqrt(Reo2) << std::endl;
    }
  }
};
