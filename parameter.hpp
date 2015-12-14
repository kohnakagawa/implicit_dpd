#pragma once

#include <string>
#include <limits>
#include <map>
#include <fstream>
#include <cassert>
#include "../../src/particle_simulator.hpp"

class Parameter {
  std::string cdir;

  void ReadTagValues(std::ifstream& fin, 
		     std::map<std::string, PS::F64>& tag_val) const 
  {
    std::string line, tag;
    double value;
    while(std::getline(fin, line)){
      std::stringstream ss(line);
      ss >> tag >> value;
      tag_val[tag] = value;
    }
  }

  void MatchingTagValues(const std::map<std::string, PS::F64>& tag_val) {
#define MATCHING(val)							\
    if(tag_val.find(std::string(#val)) == tag_val.end()){		\
      std::cerr << "Unmatching occurs." << std::endl;			\
      std::cerr << "File:" << __FILE__ << " Line:" << __LINE__ << std::endl; \
      PS::Abort();							\
    }									\
    val = tag_val[#val]							\

    MATCHING(box_leng.x);
    MATCHING(box_leng.y);
    MATCHING(box_leng.z);
    MATCHING(init_prtcl_num);
    MATCHING(dt);
    for(int i = 0; i < prop_num; i++) {
      for(int j = i; j < prop_num; j++) {
	MATCHING(cf_c[i][j]);
	MATCHING(cf_r[i][j]);
      }
    }
#undef MATCHING
  }

  void CalcGammaWithHarmonicMean(const int i, const int j) {
    if(cf_r[i][j] < 0.0){
      cf_g[i][j] = cf_g[j][i] = 2.0 / ( (1.0 / cf_g[i][i]) + (1.0 / cf_g[j][j]) );
      cf_r[i][j] = cf_r[j][i] = std::sqrt(2.0 * cf_g[i][j] * Tempera);
    }
  }
 
  void CalcInterCoef() {
    for(int i = 0; i < prop_num; i++) {
      for(int j = i + 1; j < prop_num; j++) {
	cf_c[j][i] = cf_c[i][j];
	cf_r[j][i] = cf_r[i][j];
      }
    }
    
    for(int i = 0; i < prop_num; i++)
      for(int j = 0; j < prop_num; j++)
	cf_g[i][j] = 0.5 * cf_r[i][j] * cf_r[i][j] / Tempera; // 2.0 * gamma * k_B T = sigma * sigma
    
    for(int i = 0; i < prop_num; i++)
      for(int j = i + 1; j < prop_num; j++)
	CalcGammaWithHarmonicMean(i, j);
  }

public:
  static constexpr PS::F64 Tempera = 1.0;
  
  //interactions
  static constexpr int prop_num = 2;
  PS::F64 cf_c[prop_num][prop_num]; //conservative force
  PS::F64 cf_g[prop_num][prop_num]; //dissipation force
  PS::F64 cf_r[prop_num][prop_num]; //random force
  
  //region info
  PS::F64vec box_leng, ibox_leng;
  
  //
  PS::S32 init_pnum_loc = -1, init_pnum_tot = -1;
  PS::F64 dt = std::numeric_limits<PS::F64>::quiet_NaN();
  
  Parameter(const std::string cdir_) {
    cdir = cdir_;
  }
  ~Parameter() {}

  void Initialize() {
    for(int i = 0; i < prop_num; i++) {
      for(int j = 0; j < prop_num; j++) {
	cf_c[i][j] = std::numeric_limits<PS::F64>::quiet_NaN();
	cf_g[i][j] = std::numeric_limits<PS::F64>::quiet_NaN();
	cf_r[i][j] = std::numeric_limits<PS::F64>::quiet_NaN();
	box_leng.x = box_leng.y = box_leng.z = std::numeric_limits<PS::F64>::quiet_NaN();
	ibox_leng.x = ibox_leng.y = ibox_leng.z = std::numeric_limits<PS::F64>::quiet_NaN();
      }
    }
    
  }
  
  void LoadParam() {
    const std::string fmane = cdir + "/param.txt";
    std::ifstream fin(fname.c_str());
    if(!fin) {
      std::cerr << cdir.c_str() << " dose not exist.\n";
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;
      PS::Abort();
    }
    std::map<std::string, PS::F64> tag_vals;
    ReadTagValues(fin, tag_vals);
    MatchingTagValues(tag_vals);
    CalcInterCoef();
  }

  void LoadCheck() const {
    assert(std::isfinite(box_leng.x) );
    assert(std::isfinite(box_leng.y) );
    assert(std::isfinite(box_leng.z) );

    assert(std::isfinite(ibox_leng.x) );
    assert(std::isfinite(ibox_leng.y) );
    assert(std::isfinite(ibox_leng.z) );
    
    assert(init_prtcl_num > 0);

    assert(std::isfinite(dt) );
    
    for(int i = 0; i < prop_num; i++) {
      for(int j = 0; j < prop_num; j++) {
	assert(cf_c[i][j] >= 0.0);
	assert(cf_r[i][j] >= 0.0);
	assert(cf_g[i][j] >= 0.0);

	assert(std::isfinite(cf_c[i][j]) );
	assert(std::isfinite(cf_r[i][j]) );
	assert(std::isfinite(cf_g[i][j]) );
      }
    }
  }

  void DumpAllParam() const {
    const std::string fname = cdir + "/all_param.txt";
    std::ofstream fout(fname.c_str());
#define DUMPTAGANDVAL(val)\
    fout << #val << " = " << val << std::endl \
      
    DUMPTAGANDVAL(prop_num);
    // for(int i = 0; i < prop_num; i++) {
    //   for(int j = 0; j < prop_num; j++) {
    // 	DUMPTAGANDVAL(cf_c[i][j]);
    // 	DUMPTAGANDVAL(cf_c[i][j]);
    //   }
    // }
    DUMPTAGANDVAL(box_leng.x);
    DUMPTAGANDVAL(box_leng.y);
    DUMPTAGANDVAL(box_leng.z);

    DUMPTAGANDVAL(init_prtcl_num);
    DUMPTAGANDVAL(dt);

#undef DUMPTAGANDVAL
  }
  
};
