#include <iostream>
#include <cstdio>
#include <fstream>
#include <vector>
#include <limits>
#include <map>
#include <ctime>
#include "f_calculator.hpp"
#include "parameter.hpp"

class ConfigMaker {
  std::vector<FPDPD> prtcls;

  std::string cdir, mode;
  std::ifstream fin;
  
  PS::F64 Tempera = std::numeric_limits<PS::F64>::quiet_NaN();
  PS::F64vec box_leng;
  PS::S32 amp_num = -1;
  PS::F64 lip_len = std::numeric_limits<PS::F64>::quiet_NaN();
  
  PS::F64 sph_rad = std::numeric_limits<PS::F64>::quiet_NaN();
  PS::F64 cyl_l = std::numeric_limits<PS::F64>::quiet_NaN();
  PS::F64 cyl_r = std::numeric_limits<PS::F64>::quiet_NaN();
  
  PS::F64 NormalRand(const PS::F64 mean, const PS::F64 sd) const {
    return mean + sd * std::sqrt( -2.0 * std::log(genrand_res53()) ) * std::cos(2.0 * M_PI * genrand_res53() );
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
      prtcls[i].vel = PS::F64vec(NormalRand(0.0, std::sqrt(Tempera) ),
				 NormalRand(0.0, std::sqrt(Tempera) ),
				 NormalRand(0.0, std::sqrt(Tempera) ) );
    }
    RemoveCMDrift();
  }

  void ApplyPBC(PS::F64vec& pos) {
    for(int i = 0; i < 3; i++)
      pos[i] -= std::floor(pos[i] / box_leng[i]) * box_leng[i];
  }
  
  void SetAmphilPartPos(const PS::F64vec& base, const PS::F64vec& nv, int& idx) {
    for(int unit = 0; unit < Parameter::all_unit; unit++) {
      PS::F64vec pos = base + Parameter::bond_leng * unit * nv;
      ApplyPBC(pos);
      prtcls[i].pos = pos;
      prtcls[i].id = idx;

      idx++;
    }
  }

  void MakeFlatSheet(PS::F64vec& len, PS::S32 axis) {
    assert(axis >= 0 && axis < 3);
    
    const PS::F64 eps = 1.0e-5;
    size_t idx = 0;
    PS::F64 nv(0.0, 0.0, 0.0);
    bool flap = true;
    while(idx < prtcls.size() ) {
      PS::F64vec base(len.x * genrand_res53(),
		      len.y * genrand_res53(),
		      len.z * genrand_res53());
      const PS::F64 sign = flap ? 1.0 : -1.0;
      up[axis] = len[axis] + Parameter::all_unit_n * Parameter::bond_leng + sign * eps;
      nv[axis] = -1.0 * sign;
      SetAmphilPartPos(base, nv, idx);
      flap ^= true;
    }
  }
  
  bool MakeSphLine(const PS::F64 the, const PS::F64vec& center, const PS::F64 rad,
		   const PS::F64 sign, PS::S32& idx)
  {
    const PS::F64 q_thick = 0.5 * Parameter::all_unit * Parameter::bond_leng;
    const PS::F64 prj_rad = rad * std::sin(the);

    PS::F64 d_phi = lip_len / prj_rad;
    const PS::S32 phi_elem = static_cast<int>(2.0 * M_PI / d_phi);
    d_phi = 2.0 * M_PI / phi_elem;

    for(PS::S32 i = 0; i < phi_elem; i++) {
      const PS::F64 phi = d_phi * i;
      const PS::F64vec nv(sign * std::sin(the) * std::cos(phi), 
			  sign * std::sin(the) * std::sin(phi),
			  sign * std::cos(the));
      const PS::F64vec cent2arc = sign * n_v * (rad - sign * q_thick) + center;
      SetAmphilPartPos(cent2arc, nv, idx);
      if(idx >= prtcls.size()) return false;
    }

    return true;
  }

  void MakeLineForEachTheta(const PS::S32 elem, const PS::F64 d_the, const PS::F64 offset, const PS::F64 rad, const PS::F64 sign, PS::S32& idx)
  {
    const PS::F64vec cent = 0.5 * box_leng;
    for(int i = 0; i < elem; i++) {
      const PS::F64 the = d_the * (i + offset);
      const bool flag = MakeSphLine(the, cent, rad, sign, idx);
      if(!flag) return;
    }
  }

  void MakeSphSheet() {
    const PS::F64 q_thick = 0.5 * Parameter::all_unit * Parameter::bond_leng;
    const PS::F64 out_rad = sph_rad + q_thick, in_rad  = sph_rad - q_thick;
    PS::F64 d_the_out = lip_len / out_rad, d_the_in = lip_len / in_rad;
    const PS::S32 out_the_elem = static_cast<PS::S32>(M_PI / d_the_out), in_the_elem  = static_cast<PS::S32>(M_PI / d_the_in);
    d_the_out = M_PI / out_the_elem;
    d_the_in  = M_PI / in_the_elem;
    
    PS::S32 idx = 0;
    MakeLineForEachTheta(out_the_elem, d_the_out, 0.0, out_rad, -1.0); //out
    MakeLineForEachTheta(in_the_elem, d_the_in, 0.0, in_rad, 1.0); //in
    
    const PS::S32 residue = prtcls.size() - idx;
    if(residue > 0)
      MakeLineForEachTheta(in_the_elem, d_the_in, 0.5, in_rad, -1.0);
  }

  void MakeLineForEachAxis(const PS::S32 elem, const PS::F64 d_the, const PS::F64 offset, const PS::F64 rad, const PS::F64 sign, const PS::S32 axis, PS::S32& idx)
  {
    const PS::F64 h_pi = std::acos(0.0);
    
    PS::F64vec cent = 0.5 * box_leng;
    for(int i = 0; i < elem; i++) {
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

    assert(cyl_l >= 0.0 && cyl_l < box_leng[axis]);
    assert(cyl_r >= 0.0 && (2.0 * cyl_r + Parameter::all_unit * Parameter::bond_leng) < box_leng[axis]);
    
    PS::S32 idx = 0;
    MakeLineForEachAxis(out_the_elem, d_the_out, 0.0, out_rad, -1.0); //out
    MakeLineForEachAxis(in_the_elem, d_the_in, 0.0, in_rad, 1.0);  //in
    
    const PS::S32 residue = prtcls.size() - idx;
    if(residue > 0)
      MakeLineForEachAxis(in_the_elem, d_the_in, 0.5, out_rad, -1.0);
  }

  void GenConfig() {
    switch(mode) {
    case "Plane":
      MakeFlatSheet(box_leng, 1);
      break;
    case "Sphere":
      MakeSphSheet();
      break;
    case "Cylind"
      MakeCylindSheet(2);
      break;
    default:
      if(!mode.empty() )
	std::cerr << mode << ": Unknown mode\n";
      else
	std::cerr << "Mode is not specified\n.";
      std::cerr << __FILE__ << " " __LINE__ << std::endl;
      std::exit(1);
      break;
    }
  }

  void MatchingTagValues(std::map<std::string, std::vector<std::string> >& tag_val)
  {
    Parameter::Matching(&Tempera, "Temperature", 1);
    Parameter::Matching(&(box_leng[0]), "box_leng", 3);
    Parameter::Matching(&amp_num, "amp_num", 1);
    Parameter::Matching(&mode, "mode", 1);
    if(mode == "Sphere")
      Parameter::Matching(&sph_rad, "sph_rad", 1);
    if(mode == "Cylind") {
      Parameter::Matching(&cyl_l, "cyl_l", 1);
      Parameter::Matching(&cyl_r, "cyl_r", 1);
    }
      
    Parameter::Matching(&lip_len, "lip_len", 1);
  }

public:
  ConfigMaker(const std::string cdir_) {
    cdir = cdir_;
    box_leng.x = box_leng.y = box_leng.z = std::numeric_limits<PS::F64>::quiet_NaN();
  }
  ~ConfigMaker() {
    
  }
  
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
      std::cerr << __FILE__ << " " __LINE__ std::endl;
      std::exit(1);
    }
  }

  void GenParticles() {
    prtcls.size(Parameter::all_unit * amp_num);
    GenConfig();
    GenVeloc();
  }
  
  void DumpParticleConfig() {
    std::string fname = cdir + "/init_config.txt";
    FILE* fout = fopen(fname.c_str(), "w");
    for(size_t i = 0; i < prtcls.size(); i++)
      prtcls[i].writeAscii(fout);
    fclose(fout);
  }
  
};

int main(int argc, char* argv[]) {
  if(argc != 2) {
    std::cerr << "argv[1] is target directory name.\n";
    std::exit(1);
  }
  
  ConfigMaker cmaker(std::string(argv[1]));
  cmaker.Initialize();
  cmaker.LoadParam();
  cmaker.GenParticles();
  cmaker.DumpParticleConfig();
  
  std::cout << "Initial configuration is generated at " << argv[1] << std::endl;
}
