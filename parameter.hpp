#pragma once

#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <map>
#include <cassert>
#include "io_util.hpp"

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

  void FillUpperTri(std::vector<PS::F64>& buf, PS::F64 array[], const PS::S32 dim) {
    PS::S32 cnt = 0;
    for(PS::S32 i = 0; i < dim; i++) {
      for(PS::S32 j = i; j < dim; j++) {
	array[i * dim + j] = buf[cnt];
	cnt++;
      }
    }
  }

  void MatchingTagValues(std::map<std::string, std::vector<std::string> >& tag_val) 
  {
    Matching(&(box_leng[0]) , std::string("box_leng"), tag_val, 3);
    ibox_leng.x = 1.0 / box_leng.x; ibox_leng.y = 1.0 / box_leng.y; ibox_leng.z = 1.0 / box_leng.z;
    Matching(&init_amp_num, std::string("init_amp_num"), tag_val, 1);
    Matching(&dt, std::string("dt"), tag_val, 1);
    
    std::vector<PS::F64> cf_buf(prop_num * (prop_num + 1) / 2, 0.0);
    Matching(&(cf_buf[0]), std::string("cf_r"), tag_val, prop_num * (prop_num + 1) / 2);
    FillUpperTri(cf_buf, &(cf_r[0][0]), prop_num);

    Matching(&cf_s, std::string("cf_s"), tag_val, 1);
    Matching(&cf_b, std::string("cf_b"), tag_val, 1);
    
    Matching(&chi, std::string("chi"), tag_val, 1);
    Matching(&kappa, std::string("kappa"), tag_val, 1);
    Matching(&rho_co, std::string("rho_co"), tag_val, 1);

    Matching(&all_time, std::string("all_time"), tag_val, 1);
    Matching(&step_mic, std::string("step_mic"), tag_val, 1);
    Matching(&step_mac, std::string("step_mac"), tag_val, 1);

#ifdef CHEM_MODE
    Matching(&p_thresld, std::string("p_thresld"), tag_val, 1);
    Matching(&eps, std::string("eps"), tag_val, 1);
    Matching(&max_amp_num, std::string("max_amp_num"), tag_val, 1);
#endif
  }

  void CalcGammaWithHarmonicMean(const PS::S32 i, const PS::S32 j) {
    if(cf_r[i][j] < 0.0){
      cf_g[i][j] = cf_g[j][i] = 2.0 / ( (1.0 / cf_g[i][i]) + (1.0 / cf_g[j][j]) );
      cf_r[i][j] = cf_r[j][i] = std::sqrt(2.0 * cf_g[i][j] * Tempera);
    }
  }
  
  void CalcInterCoef() {
    cf_c[Hyphob][Hyphob] = -2.0 * (kappa + 3.0) / rho_co;
    cf_c[Hyphil][Hyphil] = 0.1;
    cf_c[Hyphil][Hyphob] = cf_c[Hyphob][Hyphil] = chi / (rho_co) + 0.5 * (cf_c[Hyphil][Hyphil] + cf_c[Hyphob][Hyphob]);
    
    for(PS::S32 i = 0; i < prop_num; i++)
      for(PS::S32 j = i + 1; j < prop_num; j++)
	cf_c[j][i] = cf_c[i][j];

    cf_m[Hyphob][Hyphob][Hyphob] = 1.5 * (kappa + 2.0) / (rho_co * rho_co);
    for(PS::S32 i = 0; i < prop_num; i++)
      for(PS::S32 j= 0; j < prop_num; j++)
	for(PS::S32 k = 0; k < prop_num; k++)
	  cf_m[i][j][k] = cf_m[Hyphob][Hyphob][Hyphob];
    cf_m[Hyphil][Hyphil][Hyphil] = 0.0;
    
    for(PS::S32 i = 0; i < prop_num; i++)
      for(PS::S32 j = i + 1; j < prop_num; j++)
	cf_r[j][i] = cf_r[i][j];
    
    for(PS::S32 i = 0; i < prop_num; i++)
      for(PS::S32 j = 0; j < prop_num; j++)
	cf_g[i][j] = 0.5 * cf_r[i][j] * cf_r[i][j] / Tempera; // 2.0 * gamma * k_B T = sigma * sigma
    
    //NOTE: if cf_r[i][j] < 0.0, cf_g[i][j] and cf_r[i][j] are calculated from the harmonic mean rule.
    for(PS::S32 i = 0; i < prop_num; i++)
      for(PS::S32 j = i + 1; j < prop_num; j++)
	CalcGammaWithHarmonicMean(i, j);
    
    //NOTE: multiplied by normalize coef.
    for(PS::S32 i = 0; i < prop_num; i++) 
      for(PS::S32 j = 0; j < prop_num; j++)
	cf_r[i][j] /= std::sqrt(dt);

    const PS::F64 arc5 = arc * arc * arc * arc * arc;
    for(PS::S32 i = 0; i < prop_num; i++)
      for(PS::S32 j = 0; j < prop_num; j++)
	cf_c[i][j] *= (Reo * Reo * Reo) / (all_unit * all_unit) * (15.0 / (2.0 * M_PI * (2.0 * arc5 * arc - 3.0 * arc5 * rc + 3.0 * arc * rc2 * rc2 * rc - 2.0 * rc2 * rc2 * rc2))) * 6.0;
											 
    for(PS::S32 i = 0; i < prop_num; i++)
      for(PS::S32 j = 0; j < prop_num; j++)
	for(PS::S32 k = 0; k < prop_num; k++)
	  cf_m[i][j][k] *= ((Reo * Reo * Reo) / (all_unit * all_unit)) * (2.0 / 3.0) * (15.0 / (M_PI * rc2 * rc2)) * (15.0 / (2.0 * M_PI * rc2 * rc2 * rc)) * (Reo * Reo * Reo / all_unit);
  }

public:
  static constexpr PS::F64 Tempera	= 1.0;
  static constexpr PS::U32 head_unit	= 1;
  static constexpr PS::U32 tail_unit	= 3;
  static constexpr PS::U32 all_unit	= head_unit + tail_unit;
  static constexpr PS::F64 bond_leng	= 0.5;
  static constexpr PS::F64 ibond	= (bond_leng != 0.0) ? 1.0 / bond_leng : 0.0;
  static constexpr PS::F64 search_rad   = 1.7;
  static constexpr PS::F64 arc		= 0.9;
  static constexpr PS::F64 Reo		= 2.0;
  static constexpr PS::F64 rc           = 1.0;
  static constexpr PS::F64 rc2          = rc * rc;
  static constexpr PS::F64 irc          = 1.0 / rc;
  static constexpr PS::U32 decom_freq   = 16;

  static constexpr char atom_type[21] = {
    'O', 'N', 'C', 'S', 'P', 'Z', 'X', 'O', 'N', 'C', 'S', 'P', 'Z', 'X', 'O', 'N', 'C', 'S', 'P', 'Z', 'X'
  };
  
  enum {
    Hyphil = 0,
    Hyphob,
    
    prop_num,
  };
  
  //interactions
  static PS::F64 cf_c[prop_num][prop_num];
  static PS::F64 cf_g[prop_num][prop_num];
  static PS::F64 cf_r[prop_num][prop_num];
  static PS::F64 cf_m[prop_num][prop_num][prop_num];
  static PS::F64 cf_s;
  static PS::F64 cf_b;
  
  //region info
  static PS::F64vec box_leng, ibox_leng;
  
  PS::U32 init_amp_num = 0xffffffff, amp_num = 0xffffffff;
  PS::F64 dt = std::numeric_limits<PS::F64>::signaling_NaN();

  //macroscopic val
  PS::F64 chi = std::numeric_limits<PS::F64>::signaling_NaN();
  PS::F64 kappa = std::numeric_limits<PS::F64>::signaling_NaN();
  PS::F64 rho_co = std::numeric_limits<PS::F64>::signaling_NaN();

#ifdef CHEM_MODE
  //for chemical reaction
  PS::F64 p_thresld   = std::numeric_limits<PS::F64>::signaling_NaN();
  PS::F64 eps         = std::numeric_limits<PS::F64>::signaling_NaN();
  PS::U32 max_amp_num = 0xffffffff; //When amp_num >= max_amp_num, we stop the simulation.
#endif

#ifdef CHEM_MODE
  //for chemical reaction
  PS::F64 p_thresld = std::numeric_limits<PS::F64>::quiet_NaN();
  PS::F64 eps = std::numeric_limits<PS::F64>::quiet_NaN();
  PS::S32 max_amp_num = -1; //When amp_num >= max_amp_num, we stop the simulation.
#endif
  
  //for prng
  static PS::U32 time;
  static PS::U32 all_time, step_mic, step_mac;
  
  explicit Parameter(const std::string cdir_) {
    cdir = cdir_;
    static_assert(all_unit >= 3, "all_unit >= 3."); 
  }
  ~Parameter() {}

  template<bool ibond_is_finite>
  static inline PS::F64 cf_spring(const PS::F64 inv_dr) {
    static_assert((ibond_is_finite == true) or (ibond_is_finite == false), "bond_is_finte should be true or false");
    return;
  }

  void Initialize() {
    for(PS::S32 i = 0; i < prop_num; i++) {
      for(PS::S32 j = 0; j < prop_num; j++) {
	cf_c[i][j] = std::numeric_limits<PS::F64>::signaling_NaN();
	cf_g[i][j] = std::numeric_limits<PS::F64>::signaling_NaN();
	cf_r[i][j] = std::numeric_limits<PS::F64>::signaling_NaN();
      }
    }

    for(PS::S32 i = 0; i < prop_num; i++)
      for(PS::S32 j = 0; j < prop_num; j++)
	for(PS::S32 k = 0; k < prop_num; k++)
	  cf_m[i][j][k] = std::numeric_limits<PS::F64>::signaling_NaN();

    cf_s = cf_b = std::numeric_limits<PS::F64>::signaling_NaN();

    box_leng.x	= box_leng.y = box_leng.z = std::numeric_limits<PS::F64>::signaling_NaN();
    ibox_leng.x = ibox_leng.y = ibox_leng.z = std::numeric_limits<PS::F64>::signaling_NaN();
    
    time = 0xffffffff;
    all_time = step_mic = step_mac = 0xffffffff;
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
		       const PS::S32 num)
  {
    MatchingError(tag, tag_val, num);
    for(PS::S32 i = 0; i < num; i++)
      val[i] = std::stoi(tag_val[tag].at(i));
  }

  static void Matching(PS::U32* val,
		       std::string tag,
		       std::map<std::string, std::vector<std::string> >& tag_val,
		       const PS::U32 num)
  {
    MatchingError(tag, tag_val, num);
    for(PS::U32 i = 0; i < num; i++)
      val[i] = std::stoi(tag_val[tag].at(i));
  }

  static void Matching(PS::F64* val,
		       std::string tag,
		       std::map<std::string, std::vector<std::string> >& tag_val,
		       const PS::S32 num)
  {
    MatchingError(tag, tag_val, num);    
    for(PS::S32 i = 0; i < num; i++)
      val[i] = std::stof(tag_val[tag].at(i));
  }

  static void Matching(std::string* val,
		       std::string tag,
		       std::map<std::string, std::vector<std::string> >& tag_val,
		       const PS::S32 num)
  {
    MatchingError(tag, tag_val, num);
    for(PS::S32 i = 0; i < num; i++)
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
    const std::string fname = cdir + "/run_param.txt";
    std::ifstream fin(fname.c_str());
    if(!fin) {
      std::cerr << fname.c_str() << " dose not exist.\n";
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;
      PS::Abort();
    }
    std::map<std::string, std::vector<std::string> > tag_vals;
    ReadTagValues(fin, tag_vals);
    MatchingTagValues(tag_vals);
    CalcInterCoef();
    amp_num = init_amp_num;
  }
  
  template<class Tpsys>
  void RemoveCMDrift(Tpsys& sys) const {
    PS::F64vec cmvel(0.0, 0.0, 0.0);
    const PS::S32 num = sys.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < num; i++)
      cmvel += sys[i].vel;
    cmvel /= num;
    for(PS::S32 i = 0; i < num; i++)
      sys[i].vel -= cmvel;
  }
  
  template<class Tpsys>
  PS::U32 LoadParticleConfig(Tpsys& sys) const {
    const std::string fname = cdir + "/init_config.xyz";
    FILE* fp = io_util::xfopen(fname.c_str(), "r");
    PS::U32 line_num = 0;
    PS::U32 cur_time = 0;
    io_util::ReadXYZForm(sys, line_num, cur_time, fp);
    if(line_num / all_unit != static_cast<PS::U32>(init_amp_num) ) {
      std::cerr << line_num / all_unit << " " << init_amp_num << std::endl;
      std::cerr << "# of lines is not equal to the run input parameter information.\n";
      PS::Abort();
    }
    fclose(fp);
    RemoveCMDrift(sys);

    return cur_time;
  }

  template<class Tpsys>
  void CheckParticleConfigIsValid(const Tpsys& sys) const {
    const PS::U32 num = sys.getNumberOfParticleLocal();
    PS::F64 kin_temp = 0.0;
    for(PS::U32 i = 0; i < num; i++) {
      kin_temp += sys[i].vel * sys[i].vel;
      
      assert(sys[i].id >= 0 && sys[i].id < num);
      assert(sys[i].prop >= 0 && sys[i].prop < prop_num);
      assert(sys[i].amp_id >= 0 && sys[i].amp_id < init_amp_num);
      assert(sys[i].unit >= 0 && sys[i].unit < all_unit);
      
      for(PS::U32 j = 0; j < 3; j++) {
	if(!(sys[i].pos[j] <= box_leng[j] && sys[i].pos[j] >= 0.0)) {
	  std::cerr << "There is a particle in outside range.\n";
	  std::cerr << __FILE__ << " " << __LINE__ << std::endl;
	  PS::Abort();
	}
      }
    }

    kin_temp /= 3.0 * num;
    if(fabs(kin_temp - Tempera) >= 1.0e-1) {
      std::cerr << "Please check kinetic temperature.\n";
      std::cerr << "Kinetic temperature is " << kin_temp << "K_BT.\n";
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;
      PS::Abort();
    }
  }

  void CheckLoaded() const {
    for(PS::S32 i = 0; i < prop_num; i++) {
      for(PS::S32 j = 0; j < prop_num; j++) {
	assert(std::isfinite(cf_c[i][j]) );
	assert(std::isfinite(cf_r[i][j]) );
	assert(std::isfinite(cf_g[i][j]) );
      }
    }

    for(PS::S32 i = 0; i < prop_num; i++)
      for(PS::S32 j = 0; j < prop_num; j++)
	for(PS::S32 k = 0; k < prop_num; k++)
	  assert(std::isfinite(cf_m[i][j][k]));

    assert(std::isfinite(cf_s));
    assert(std::isfinite(cf_b));

    assert(std::isfinite(box_leng.x) );
    assert(std::isfinite(box_leng.y) );
    assert(std::isfinite(box_leng.z) );

    assert(std::isfinite(ibox_leng.x) );
    assert(std::isfinite(ibox_leng.y) );
    assert(std::isfinite(ibox_leng.z) );

    assert(init_amp_num > 0);
    assert(amp_num > 0);

    assert(std::isfinite(dt) );
    
    assert(std::isfinite(chi) );
    assert(std::isfinite(kappa));
    assert(std::isfinite(rho_co));
    
#ifdef CHEM_MODE
    assert(std::isfinite(p_thresld));
    assert(p_thresld <= 1.0 && p_thresld >= 0.0);
    assert(std::isfinite(eps));
    assert(max_amp_num >= init_amp_num);
    assert(max_amp_num * Parameter::all_unit < std::numeric_limits<PS::U32>::max() );
    assert(max_amp_num != 0xffffffff);
#endif

    assert(time != 0xffffffff);
    assert(all_time != 0xffffffff);
    assert(step_mic != 0xffffffff);
    assert(step_mac != 0xffffffff);
  }

  void DumpAllParam(std::ostream& ost) const {
#define DUMPTAGANDVAL(val) ost << #val << " = " << val << std::endl

    DUMPTAGANDVAL(cdir);
    
    DUMPTAGANDVAL(Tempera);
    DUMPTAGANDVAL(head_unit);
    DUMPTAGANDVAL(tail_unit);
    DUMPTAGANDVAL(all_unit);
    DUMPTAGANDVAL(bond_leng);
    DUMPTAGANDVAL(ibond);
    DUMPTAGANDVAL(search_rad);
    DUMPTAGANDVAL(arc);
    DUMPTAGANDVAL(Reo);
    
    DUMPTAGANDVAL(prop_num);
    DUMPTAGANDVAL(rc);
    DUMPTAGANDVAL(rc2);
    DUMPTAGANDVAL(irc);
    
    DUMPTAGANDVAL(box_leng.x);
    DUMPTAGANDVAL(box_leng.y);
    DUMPTAGANDVAL(box_leng.z);

    DUMPTAGANDVAL(ibox_leng.x);
    DUMPTAGANDVAL(ibox_leng.y);
    DUMPTAGANDVAL(ibox_leng.z);

    DUMPTAGANDVAL(init_amp_num);
    DUMPTAGANDVAL(amp_num);
    DUMPTAGANDVAL(dt);

    DUMPTAGANDVAL(chi);
    DUMPTAGANDVAL(kappa);
    DUMPTAGANDVAL(rho_co);

    ost << "NOTE:\n";
    ost << "cf_r are multiplied by 1 / sqrt(dt).\n";
      
#define DUMPINTRPARAM(val) ost << #val << ":\n";		\
    for(PS::S32 i = 0; i < prop_num; i++) {			\
      for(PS::S32 j = 0; j < prop_num; j++)			\
	ost << val[i][j] << " ";				\
      ost << std::endl;					\
    } 

    DUMPINTRPARAM(cf_c);
    DUMPINTRPARAM(cf_r);
    DUMPINTRPARAM(cf_g);
    
    ost << "cf_m:\n";
    for(PS::S32 i = 0; i < prop_num; i++)
      for(PS::S32 j = 0; j < prop_num; j++)
	for(PS::S32 k = 0; k < prop_num; k++)
	  ost << cf_m[i][j][k] << " ";
    ost << std::endl;

    DUMPTAGANDVAL(cf_s);
    DUMPTAGANDVAL(cf_b);

#undef DUMPTAGANDVAL
#undef DUMPINTRPARAM
  }

  void DumpAllParam() const {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    if(PS::Comm::getRank() == 0) {
#endif
      const std::string fname = cdir + "/all_param.txt";
      std::ofstream fout(fname.c_str());
      DumpAllParam(fout);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    }
#endif
  }

  //for share data
  void ShareDataWithOtherProc() {
    const PS::S32 num_two_body = Parameter::prop_num * Parameter::prop_num;
    const PS::S32 num_thre_body = Parameter::prop_num * Parameter::prop_num * Parameter::prop_num;
    
    // static member
    PS::Comm::broadcast(&(cf_c[0][0]), num_two_body, 0);
    PS::Comm::broadcast(&(cf_g[0][0]), num_two_body, 0);
    PS::Comm::broadcast(&(cf_r[0][0]), num_two_body, 0);
    PS::Comm::broadcast(&(cf_m[0][0][0]), num_thre_body, 0);
    PS::Comm::broadcast(&cf_s, Parameter::prop_num, 0);
    PS::Comm::broadcast(&cf_b, Parameter::prop_num, 0);
    
    PS::Comm::broadcast(&box_leng, 1, 0);
    PS::Comm::broadcast(&ibox_leng, 1, 0);
    
    PS::Comm::broadcast(&time, 1, 0);
    PS::Comm::broadcast(&all_time, 1, 0);
    PS::Comm::broadcast(&step_mic, 1, 0);
    PS::Comm::broadcast(&step_mac, 1, 0);
    
    // non static member
    PS::Comm::broadcast(&init_amp_num, 1, 0);
    PS::Comm::broadcast(&amp_num, 1, 0);
    PS::Comm::broadcast(&dt, 1, 0);
    PS::Comm::broadcast(&chi, 1, 0);
    PS::Comm::broadcast(&kappa, 1, 0);
    PS::Comm::broadcast(&rho_co, 1, 0);

#ifdef CHEM_MODE
    PS::Comm::broadcast(&p_thresld, 1, 0);
    PS::Comm::broadcast(&eps, 1, 0);
    PS::Comm::broadcast(&max_amp_num, 1, 0);
#endif
  }
  
  //for debug
  void DebugDumpAllParam() const {
    std::stringstream ss;
    const PS::S32 rank = PS::Comm::getRank();
    ss << cdir << "/debug_dump_rank" << rank << ".txt";

    std::ofstream fout(ss.str().c_str());
    DumpAllParam(fout);
  }
};

template<>
inline PS::F64 Parameter::cf_spring<true>(const PS::F64 inv_dr) {
  return Parameter::cf_s * (inv_dr - Parameter::ibond);
}

template<>
inline PS::F64 Parameter::cf_spring<false>(const PS::F64 inv_dr) {
  return -Parameter::cf_s;
}
