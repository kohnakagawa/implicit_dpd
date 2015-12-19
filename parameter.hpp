#pragma once

#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <map>
#include <cassert>

class Parameter {
  std::string cdir;

  static void DeleteHeadSpace(std::string &buf) {
    while(buf.find_first_of(" 　\t") == 0) {
      buf.erase(buf.begin());
      if(buf.empty()) break;
    }
  }

  static void DeleteTailSpace(std::string &buf) {
    size_t tail = buf.size() - 1;
    while(buf.find_last_of(" 　\t") == tail) {
      buf.erase(buf.end() - 1);
      if(buf.empty()) break;
      tail = buf.size() - 1;
    }
  }

  void FillUpperTri(std::vector<PS::F64>& buf, PS::F64 array[], const int dim) {
    int cnt = 0;
    for(int i = 0; i < dim; i++) {
      for(int j = i; j < dim; j++) {
	array[i * dim + j] = buf[cnt];
	cnt++;
      }
    }
  }

  void MatchingTagValues(std::map<std::string, std::vector<std::string> >& tag_val) 
  {
    Matching(&(box_leng[0]) , std::string("box_leng"), tag_val, 3);
    Matching(&init_prtcl_num, std::string("init_prtcl_num"), tag_val, 1);
    Matching(&dt, std::string("dt"), tag_val, 1);
    
    std::vector<PS::F64> cf_buf(prop_num * (prop_num + 1) / 2, 0.0);
    Matching(&(cf_buf[0]), std::string("cf_c"), tag_val, prop_num * (prop_num + 1) / 2);
    FillUpperTri(cf_buf, &(cf_c[0][0]), prop_num);
    Matching(&(cf_buf[0]), std::string("cf_r"), tag_val, prop_num * (prop_num + 1) / 2);
    FillUpperTri(cf_buf, &(cf_r[0][0]), prop_num);

    Matching(&rc, std::string("rc"), tag_val, 1);
    rc2 = rc * rc;
    irc = 1.0 / rc;

    Matching(&cf_s, std::string("cf_s"), tag_val, 1);
    Matching(&cf_b, std::string("cf_b"), tag_val, 1);
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

    for(int i = 0; i < prop_num; i++) 
      for(int j = 0; j < prop_num; j++)
	cf_r[i][j] /= std::sqrt(dt);
  }

public:
  static constexpr PS::F64 Tempera	= 1.0;
  static constexpr PS::S32 head_unit	= 1;
  static constexpr PS::S32 tail_unit	= 3;
  static constexpr PS::S32 all_unit	= head_unit + tail_unit;
  static constexpr PS::F64 bond_leng	= 0.5;
  static constexpr PS::F64 ibond	= 1.0 / bond_leng;
  static constexpr PS::S32 prop_num = 2;
  
  //interactions
  static PS::F64 cf_c[prop_num][prop_num];
  static PS::F64 cf_g[prop_num][prop_num];
  static PS::F64 cf_r[prop_num][prop_num];
  static PS::F64 cf_s;
  static PS::F64 cf_b;

  //cutoff length
  static PS::F64 rc, rc2, irc;
  
  //region info
  static PS::F64vec box_leng, ibox_leng;
  
  PS::S32 init_prtcl_num = -1;
  PS::F64 dt = std::numeric_limits<PS::F64>::quiet_NaN();

  //for prng
  static PS::U32 time;
  
  Parameter(const std::string& cdir_) { 
    cdir = cdir_;
    static_assert(all_unit >= 3, "all_unit >= 3."); 
    cf_s = cf_b = rc = rc2 = irc = std::numeric_limits<PS::F64>::quiet_NaN();
    time = 0xffffffff;
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

  static void MatchingError(std::string tag, 
			    std::map<std::string, std::vector<std::string> >& tag_val,
			    const size_t num) {
    if(tag_val.find(tag) == tag_val.end()){				
      std::cerr << "Unmatching occurs." << std::endl;			
      std::cerr << "File:" << __FILE__ << " Line:" << __LINE__ << std::endl;
      PS::Abort();							
    }									
    assert(tag_val[tag].size() == num);					
  }

  static void Matching(PS::S32* val,
		       std::string tag,
		       std::map<std::string, std::vector<std::string> >& tag_val,
		       const int num)
  {
    MatchingError(tag, tag_val, num);
    for(int i = 0; i < num; i++)
      val[i] = std::stoi(tag_val[tag].at(i));
  }

  static void Matching(PS::F64* val,
		       std::string tag,
		       std::map<std::string, std::vector<std::string> >& tag_val,
		       const int num)
  {
    MatchingError(tag, tag_val, num);    
    for(int i = 0; i < num; i++)
      val[i] = std::stof(tag_val[tag].at(i));
  }

  static void Matching(std::string* val,
		       std::string tag,
		       std::map<std::string, std::vector<std::string> >& tag_val,
		       const int num)
  {
    MatchingError(tag, tag_val, num);
    for(int i = 0; i < num; i++)
      val[i] = tag_val[tag].at(i);
  }

  static void ReadTagValues(std::ifstream& fin, 
			    std::map<std::string, std::vector<std::string> >& tag_val)
  {
    std::string line, tag;
    std::string::size_type comment_start = 0;
    while(std::getline(fin, line) ) {
      if( (comment_start = line.find(';')) != std::string::size_type(-1) )
	line = line.substr(0, comment_start);

      if(line.empty() ) continue;

      DeleteHeadSpace(line);
      DeleteTailSpace(line);
      
      std::stringstream ss(line);
      ss >> tag;
      std::vector<std::string> values;
      while(!ss.eof() ) {
	std::string buf;
	ss >> buf;
	values.push_back(buf);
      }
      tag_val[tag] = values;
    }
  }
  
  void LoadParam() {
    const std::string fname = cdir + "/param.txt";
    std::ifstream fin(fname.c_str());
    if(!fin) {
      std::cerr << cdir.c_str() << " dose not exist.\n";
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;
      PS::Abort();
    }
    std::map<std::string, std::vector<std::string> > tag_vals;
    ReadTagValues(fin, tag_vals);
    MatchingTagValues(tag_vals);
    CalcInterCoef();
  }

  template<class Tpsys>
  bool LoadParticleConfig(Tpsys& sys) const {
    const std::string fname = cdir + "/init_config.txt";
    if(fopen(fname.c_str(), "r") != NULL) {
      sys.readParticleAscii(fname.c_str());
      return true;
    } else {
      std::cerr << "init_config.txt does not exist.\n";
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;
      return false;
    }
  }

  template<class Tpsys>
  void CheckParticleConfigIsValid(const Tpsys& sys) const {
    const PS::S32 num = sys.getNumberOfParticleLocal();
    PS::F64 kin_temp = 0.0;
    for(PS::S32 i = 0; i < num; i++) {
      kin_temp += sys[i].vel * sys[i].vel;
      for(PS::S32 j = 0; j < 3; j++) {
	if(!(sys[i].pos[j] <= box_leng[j] && sys[i].pos[j] >= 0.0)) {
	  std::cerr << "There is a particle in outside range.\n";
	  std::cerr << __FILE__ << " " << __LINE__ << std::endl;
	  PS::Abort();
	}
      }
    }

    kin_temp /= 3.0 * num;
    if(fabs(kin_temp - Tempera) >= 1.0e-2) {
      std::cerr << "Please check kinetic temperature.\n";
      std::cerr << "Kinetic temperature is " << kin_temp << "K_BT.\n";
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;
      PS::Abort();
    }
  }

  void CheckLoaded() const {
    assert(std::isfinite(box_leng.x) );
    assert(std::isfinite(box_leng.y) );
    assert(std::isfinite(box_leng.z) );

    assert(std::isfinite(ibox_leng.x) );
    assert(std::isfinite(ibox_leng.y) );
    assert(std::isfinite(ibox_leng.z) );

    assert(std::isfinite(rc) );
    assert(std::isfinite(rc2) );
    assert(std::isfinite(irc) );
    
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
    
#define DUMPTAGANDVAL(val) fout << #val << " = " << val << std::endl

    DUMPTAGANDVAL(prop_num);
    DUMPTAGANDVAL(rc);
    DUMPTAGANDVAL(rc2);
    DUMPTAGANDVAL(irc);
    
    DUMPTAGANDVAL(box_leng.x);
    DUMPTAGANDVAL(box_leng.y);
    DUMPTAGANDVAL(box_leng.z);

    DUMPTAGANDVAL(init_prtcl_num);
    DUMPTAGANDVAL(dt);

#undef DUMPTAGANDVAL    

    fout << "NOTE:\n";
    fout << "cf_r are multiplied by 1 / sqrt(dt).";
      
#define DUMPINTRPARAM(val) fout << #val << ":\n";	\
    for(int i = 0; i < prop_num; i++) {			\
      for(int j = 0; j < prop_num; j++)			\
	fout << val[i][j] << " ";			\
      fout << std::endl;				\
    } 

    DUMPINTRPARAM(cf_c);
    DUMPINTRPARAM(cf_r);
    DUMPINTRPARAM(cf_g);

#undef DUMPINTRPARAM
  }
};

PS::F64 Parameter::cf_c[Parameter::prop_num][Parameter::prop_num];
PS::F64 Parameter::cf_g[Parameter::prop_num][Parameter::prop_num];
PS::F64 Parameter::cf_r[Parameter::prop_num][Parameter::prop_num];
PS::F64 Parameter::cf_s;
PS::F64 Parameter::cf_b;

PS::F64 Parameter::rc, Parameter::rc2, Parameter::irc;
PS::F64vec Parameter::box_leng, Parameter::ibox_leng;

PS::U32 Parameter::time;
