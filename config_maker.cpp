#include <iostream>
#include <cstdio>
#include <vector>
#include <ctime>
#include "particle_simulator.hpp"
#include "io_util.hpp"
#include "parameter.hpp"
#include "f_calculator.hpp"
#include <cstring>

static_assert(Parameter::bond_leng != 0.0, "Please check Parameter::bond_leng != 0.0");

class ConfigMaker {
  std::vector<FPDPD> prtcls;

  std::string cdir, mode;
  std::ifstream fin;
  
  PS::F64 Tempera	= std::numeric_limits<PS::F64>::quiet_NaN();
  PS::F64vec box_leng;
  PS::U32 amp_num	= 0xffffffff;
  PS::F64 lip_len	= std::numeric_limits<PS::F64>::quiet_NaN();
  
  PS::F64 sph_rad	= std::numeric_limits<PS::F64>::quiet_NaN();
  PS::F64 cyl_l		= std::numeric_limits<PS::F64>::quiet_NaN();
  PS::F64 cyl_r		= std::numeric_limits<PS::F64>::quiet_NaN();
  
  PS::F64 NormalRand(const PS::F64 mean, const PS::F64 sd) const {
    return mean + sd * std::sqrt( -2.0 * std::log(PS::MT::genrand_real3()) ) * std::cos(2.0 * M_PI * PS::MT::genrand_real3() );
  }
  
  void RemoveCMDrift() {
    PS::F64vec cmvel(0.0, 0.0, 0.0);
    for(size_t i = 0; i < prtcls.size(); i++)
      cmvel += prtcls[i].vel;
    cmvel /= prtcls.size();
    for(size_t i = 0; i < prtcls.size(); i++)
      prtcls[i].vel -= cmvel;
  }

  void GenVeloc() {
    for(size_t i = 0; i < prtcls.size(); i++) {
      prtcls[i].vel =  PS::F64vec(NormalRand(0.0, std::sqrt(Tempera) ),
				  NormalRand(0.0, std::sqrt(Tempera) ),
				  NormalRand(0.0, std::sqrt(Tempera) ) );
    }
    RemoveCMDrift();
    
    for(size_t i = 0; i < prtcls.size(); i++) {
      prtcls[i].vel_buf = prtcls[i].vel;
      prtcls[i].acc = PS::F64vec(0.0, 0.0, 0.0);
    }
  }

  void ApplyPBC(PS::F64vec& pos) {
    for(PS::U32 i = 0; i < 3; i++)
      pos[i] -= std::floor(pos[i] / box_leng[i]) * box_leng[i];
  }
  
  void SetAmphilPartPos(const PS::F64vec& base, const PS::F64vec& nv, PS::U32& idx) {
    for(PS::U32 unit = 0; unit < Parameter::all_unit; unit++) {
      PS::F64vec pos = base + Parameter::bond_leng * unit * nv;
      ApplyPBC(pos);
      prtcls[idx].pos = pos;
      prtcls[idx].id = idx;
      prtcls[idx].unit = unit;
      prtcls[idx].amp_id = idx / Parameter::all_unit;
      if(unit < Parameter::head_unit)
	prtcls[idx].prop = Parameter::Hyphil;
      else
	prtcls[idx].prop = Parameter::Hyphob;
      idx++;
    }
  }

  void MakeFlatSheet(PS::F64vec& len, PS::S32 axis) {
    assert(axis >= 0 && axis < 3);
    
    const PS::F64 eps = 1.0e-5;
    PS::U32 idx = 0;
    PS::F64vec nv(0.0, 0.0, 0.0);
    bool flap = true;
    while(idx < prtcls.size() ) {
      PS::F64vec base(len.x * PS::MT::genrand_res53(),
		      len.y * PS::MT::genrand_res53(),
		      len.z * PS::MT::genrand_res53());
      const PS::F64 sign = flap ? 1.0 : -1.0;
      base[axis] = 0.5 * len[axis] + (Parameter::all_unit * Parameter::bond_leng + eps) * sign;
      nv[axis] = -1.0 * sign;
      SetAmphilPartPos(base, nv, idx);
      flap ^= true;
    }
  }
  
  bool MakeSphLine(const PS::F64 the, const PS::F64vec& center, const PS::F64 rad,
		   const PS::F64 sign, PS::U32& idx)
  {
    const PS::F64 q_thick = 0.5 * Parameter::all_unit * Parameter::bond_leng;
    const PS::F64 prj_rad = rad * std::sin(the);

    PS::F64 d_phi = lip_len / prj_rad;
    const PS::S32 phi_elem = static_cast<PS::S32>(2.0 * M_PI / d_phi);
    d_phi = 2.0 * M_PI / phi_elem;

    for(PS::S32 i = 0; i < phi_elem; i++) {
      const PS::F64 phi = d_phi * i;
      const PS::F64vec nv(sign * std::sin(the) * std::cos(phi), 
			  sign * std::sin(the) * std::sin(phi),
			  sign * std::cos(the));
      const PS::F64vec cent2arc = sign * nv * (rad - sign * q_thick) + center;
      SetAmphilPartPos(cent2arc, nv, idx);
      if(idx >= prtcls.size()) return false;
    }

    return true;
  }

  void MakeLineForEachTheta(const PS::S32 elem,
			    const PS::F64 d_the,
			    const PS::F64 offset,
			    const PS::F64 rad,
			    const PS::F64 sign,
			    PS::U32& idx)
  {
    const PS::F64vec cent = 0.5 * box_leng;
    for(PS::S32 i = 0; i < elem; i++) {
      const PS::F64 the = d_the * (i + offset);
      const bool flag = MakeSphLine(the, cent, rad, sign, idx);
      if(!flag) return;
    }
  }

  void MakeSphSheet() {
    const PS::F64 min_box_leng = box_leng.getMin();
    assert(sph_rad >= 0.0 && sph_rad <= 0.5 * min_box_leng);
    
    const PS::F64 q_thick = 0.5 * Parameter::all_unit * Parameter::bond_leng;
    const PS::F64 out_rad = sph_rad + q_thick, in_rad  = sph_rad - q_thick;
    PS::F64 d_the_out = lip_len / out_rad, d_the_in = lip_len / in_rad;
    const PS::S32 out_the_elem = static_cast<PS::S32>(M_PI / d_the_out), in_the_elem  = static_cast<PS::S32>(M_PI / d_the_in);
    d_the_out = M_PI / out_the_elem;
    d_the_in  = M_PI / in_the_elem;
    
    PS::U32 idx = 0;
    MakeLineForEachTheta(out_the_elem, d_the_out, 0.0, out_rad, -1.0, idx); //out
    MakeLineForEachTheta(in_the_elem, d_the_in, 0.0, in_rad, 1.0, idx); //in
    
    const PS::U32 residue = prtcls.size() - idx;
    if(residue > 0)
      MakeLineForEachTheta(in_the_elem, d_the_in, 0.5, in_rad, -1.0, idx);
  }

  void MakeLineForEachAxis(const PS::S32 elem,
			   const PS::F64 d_the,
			   const PS::F64 offset,
			   const PS::F64 rad,
			   const PS::F64 sign,
			   const PS::S32 axis,
			   PS::U32& idx)
  {
    const PS::F64 h_pi = std::acos(0.0);
    
    PS::F64vec cent = 0.5 * box_leng;
    for(PS::S32 i = 0; i < elem; i++) {
      cent[axis] = elem * lip_len;
      const bool flag = MakeSphLine(h_pi, cent, rad, sign, idx);
      if(!flag) return;
    }
  }

  void MakeCylindSheet(const PS::S32 axis) {
    assert(axis >= 0 && axis < 3);
    
    const PS::F64 q_thick = 0.5 * Parameter::all_unit * Parameter::bond_leng;
    const PS::F64 out_rad = cyl_r + q_thick, in_rad  = cyl_r - q_thick;
    PS::F64 d_the_out = lip_len / out_rad, d_the_in = lip_len / in_rad;
    const PS::S32 out_the_elem = static_cast<PS::S32>(M_PI / d_the_out), in_the_elem  = static_cast<PS::S32>(M_PI / d_the_in);
    d_the_out = M_PI / out_the_elem; d_the_in  = M_PI / in_the_elem;

    PS::F64vec cent = 0.5 * box_leng;
    cent[axis] = 0.0;

    assert(cyl_l >= 0.0 && cyl_l <= box_leng[axis]);
    assert(cyl_r >= 0.0 && (2.0 * cyl_r + Parameter::all_unit * Parameter::bond_leng) < box_leng[axis]);
    
    PS::U32 idx = 0;
    MakeLineForEachAxis(out_the_elem, d_the_out, 0.0, out_rad, -1.0, axis, idx); //out
    MakeLineForEachAxis(in_the_elem, d_the_in, 0.0, in_rad, 1.0, axis, idx);  //in
    
    const PS::U32 residue = prtcls.size() - idx;
    if(residue > 0)
      MakeLineForEachAxis(in_the_elem, d_the_in, 0.5, out_rad, -1.0, axis, idx);
  }

  void MakeRandomConfig(PS::F64vec& len) {
    PS::U32 idx = 0;
    while(idx < prtcls.size() ) {
      PS::F64vec base(len.x * PS::MT::genrand_res53(),
		      len.y * PS::MT::genrand_res53(),
		      len.z * PS::MT::genrand_res53());
      
      PS::F64vec nv(PS::MT::genrand_res53(),
		    PS::MT::genrand_res53(),
		    PS::MT::genrand_res53());
      
      const PS::F64 norm_nv = std::sqrt(nv * nv);
      nv /= norm_nv;
      SetAmphilPartPos(base, nv, idx);
    }
  }

#define MODE_EQ(str) strcasecmp(mode.c_str(), str) == 0
  void GenConfig() {
    if(MODE_EQ("flat")) {
      MakeFlatSheet(box_leng, 1);      
    } else if(MODE_EQ("sphere")) {
      MakeSphSheet();
    } else if(MODE_EQ("cylind")) {
      MakeCylindSheet(2);
    } else if(MODE_EQ("random")) {
      MakeRandomConfig(box_leng);
    } else {
      if(!mode.empty() )
	std::cerr << mode << ": Unknown mode\n";
      else
	std::cerr << "Mode is not specified\n.";
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;
      std::exit(1);
    }
  }

  void MatchingTagValues(std::map<std::string, std::vector<std::string> >& tag_val)
  {
    Parameter::Matching(&Tempera, std::string("Temperature"), tag_val, 1);
    Parameter::Matching(&(box_leng[0]), std::string("box_leng"), tag_val, 3);
    Parameter::Matching(&amp_num, "amp_num", tag_val, 1);
    Parameter::Matching(&mode, "mode", tag_val, 1);
    if(MODE_EQ("sphere"))
      Parameter::Matching(&sph_rad, "sph_rad", tag_val, 1);
    if(MODE_EQ("cylind")) {
      Parameter::Matching(&cyl_l, "cyl_l", tag_val, 1);
      Parameter::Matching(&cyl_r, "cyl_r", tag_val, 1);
    }
      
    Parameter::Matching(&lip_len, "lip_len", tag_val, 1);
  }
#undef MODE_EQ
  
  void InitializeParticle() {
    for(size_t i = 0; i < prtcls.size(); i++) {
      prtcls[i].id = prtcls[i].prop = prtcls[i].amp_id = prtcls[i].unit = 0xffffffff;
      prtcls[i].pos.x = prtcls[i].pos.y = prtcls[i].pos.z = std::numeric_limits<PS::F64>::quiet_NaN();
      prtcls[i].vel.x = prtcls[i].vel.y = prtcls[i].vel.z = std::numeric_limits<PS::F64>::quiet_NaN();
      prtcls[i].vel_buf.x = prtcls[i].vel_buf.y = prtcls[i].vel_buf.z = 0.0;
      prtcls[i].acc.x = prtcls[i].acc.y = prtcls[i].acc.z = 0.0;
    }
  }

  void CheckParticleConfigIsValid() {
    PS::F64 kin_temp = 0.0;
    for(size_t i = 0; i < prtcls.size(); i++) {
      kin_temp += prtcls[i].vel * prtcls[i].vel;
      for(int j = 0; j < 3; j++) {
	if(!(prtcls[i].pos[j] <= box_leng[j] && prtcls[i].pos[j] >= 0.0) ) {
	  std::cerr << "There is a particle in outside range.\n";
	  std::cerr << __FILE__ << " " << __LINE__ << std::endl;
	  PS::Abort();
	}
      }
    }
    
    kin_temp /= 3.0 * prtcls.size();
    if(fabs(kin_temp - Tempera) >= 1.0e-2) {
      std::cerr << "Please check kinetic temperature.\n";
      std::cerr << "Kinetic temperature is " << kin_temp << "K_BT.\n";
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;
      PS::Abort();
    }

    const PS::U32 all_n = Parameter::all_unit * amp_num;
    for(size_t i = 0; i < prtcls.size(); i++) {
      assert(prtcls[i].id >= 0 && prtcls[i].id < all_n);
      assert(prtcls[i].prop >= 0 && prtcls[i].prop < Parameter::prop_num);
      assert(prtcls[i].amp_id >= 0 && prtcls[i].amp_id < amp_num);
      assert(prtcls[i].unit >= 0 && prtcls[i].unit < Parameter::all_unit);
    }
  }
  
public:
  explicit ConfigMaker(const std::string& cdir_) {
    cdir = cdir_;
    box_leng.x = box_leng.y = box_leng.z = std::numeric_limits<PS::F64>::quiet_NaN();
  }
  ~ConfigMaker() {}
  
  void Initialize() {
    PS::MT::init_genrand( static_cast<unsigned long>(time(NULL)) );
  }

  void LoadParam() {
    std::string fname = cdir + "/config_param.txt";
    std::ifstream fin(fname.c_str());
    if(fin) {
      std::map<std::string, std::vector<std::string> > tag_val;
      Parameter::ReadTagValues(fin, tag_val);
      MatchingTagValues(tag_val);
    } else {
      std::cerr << "Cannot find file " << fname << std::endl;
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;
      std::exit(1);
    }
  }

  void GenParticles() {
    const PS::U32 all_n = Parameter::all_unit * amp_num;
    assert(all_n > 0 && all_n < std::numeric_limits<PS::U32>::max());
    prtcls.resize(all_n);
    InitializeParticle();
    GenConfig();
    GenVeloc();
    CheckParticleConfigIsValid();
  }
  
  void DumpParticleConfig() {
    std::string fname = cdir + "/init_config.xyz";
    FILE* fout = io_util::xfopen(fname.c_str(), "w");
    io_util::WriteXYZForm(&prtcls[0], prtcls.size(), 0, fout);
    fclose(fout);
  } 
};

int main(int argc, char* argv[]) {
  if(argc != 2) {
    std::cerr << "argv[1] is target directory name.\n";
    std::exit(1);
  }

  const std::string cdir = argv[1];
  ConfigMaker cmaker(cdir);
  cmaker.Initialize();
  cmaker.LoadParam();
  cmaker.GenParticles();
  cmaker.DumpParticleConfig();
  
  std::cout << "Initial configuration is generated at " << argv[1] << std::endl;
}
