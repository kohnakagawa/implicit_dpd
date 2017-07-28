#include <iostream>
#include <cstdio>
#include <vector>
#include <ctime>
#include "particle_simulator.hpp"
#include "io_util.hpp"
#include "parameter.hpp"
#include "f_calculator.hpp"
#include <cstring>

constexpr char Parameter::atom_type[21];

static_assert(Parameter::bond_leng != 0.0, "Please check Parameter::bond_leng != 0.0");

class ConfigMaker {
  std::vector<FPDPD> prtcls;

  std::string cdir, mode;
  std::ifstream fin;

  PS::F64 Tempera      = std::numeric_limits<PS::F64>::signaling_NaN();
  PS::F64vec box_leng;
  PS::U32 amp_num      = 0xffffffff;
  PS::U32 amp_ptcl_num = 0xffffffff;
  PS::F64 lip_len      = std::numeric_limits<PS::F64>::signaling_NaN();

  PS::F64 sph_rad      = std::numeric_limits<PS::F64>::signaling_NaN();
  PS::F64 upper_the    = std::numeric_limits<PS::F64>::signaling_NaN();

  PS::F64 cyl_l        = std::numeric_limits<PS::F64>::signaling_NaN();
  PS::F64 cyl_r        = std::numeric_limits<PS::F64>::signaling_NaN();

  PS::F64 in_out_rat   = std::numeric_limits<PS::F64>::signaling_NaN();

  PS::U32 sol_num      = 0xffffffff;
  bool with_solvent    = false;

  static PS::F64 NormalRand(const PS::F64 mean, const PS::F64 sd) {
    return mean + sd * std::sqrt( -2.0 * std::log(PS::MT::genrand_real3())) * std::cos(2.0 * M_PI * PS::MT::genrand_real3());
  }

  static PS::F64vec GenrandF64vec(const PS::F64vec& vec) {
    const PS::F64vec rndvec(vec.x * PS::MT::genrand_real3(),
                            vec.y * PS::MT::genrand_real3(),
                            vec.z * PS::MT::genrand_real3());
    return rndvec;
  }

  void RemoveCMDrift() {
    PS::F64vec cmvel(0.0, 0.0, 0.0);
    for (const auto& prtcl : prtcls)
      cmvel += prtcl.vel;
    cmvel /= prtcls.size();
    for (auto& prtcl : prtcls)
      prtcl.vel -= cmvel;
  }

  void GenVeloc() {
    for (auto& prtcl : prtcls) {
      prtcl.vel = PS::F64vec(NormalRand(0.0, std::sqrt(Tempera)),
                             NormalRand(0.0, std::sqrt(Tempera)),
                             NormalRand(0.0, std::sqrt(Tempera)));
    }
    RemoveCMDrift();

    for (auto& prtcl : prtcls) {
      prtcl.vel_buf = prtcl.vel;
      prtcl.acc = PS::F64vec(0.0, 0.0, 0.0);
    }
  }

  void ApplyPBC(PS::F64vec& pos) {
    for(PS::U32 i = 0; i < 3; i++)
      pos[i] -= std::floor(pos[i] / box_leng[i]) * box_leng[i];
  }

  void SetAmphilPartPos(const PS::F64vec& base, const PS::F64vec& nv, PS::U32& idx) {
    for (PS::U32 unit = 0; unit < Parameter::all_unit; unit++) {
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

  void SetSolventPos(const PS::F64vec& pos, const PS::U32 idx) {
    prtcls[idx].pos = pos;
    prtcls[idx].id  = idx;
    prtcls[idx].prop = Parameter::Solvent;
    prtcls[idx].unit = 0xffffffff;
    prtcls[idx].amp_id = 0xffffffff;
  }

  void MakeFlatSheet(PS::F64vec& len,
                     const PS::S32 axis,
                     bool asym,
                     const PS::U32 num,
                     const PS::U32 offset) {
    assert(axis >= 0 && axis < 3);
    PS::U32 rat = 0xffffffff;
    if (asym) {
      assert(mode == "asym_flat");
      rat = static_cast<PS::U32>(100 * in_out_rat);
    } else {
      assert(mode == "flat");
    }

    const PS::F64 eps = 1.0e-5;
    PS::U32 idx = offset;
    const PS::U32 end_idx = offset + num;
    PS::F64vec nv(0.0, 0.0, 0.0);
    bool flap = true;
    while (idx < end_idx) {
      PS::F64vec base = GenrandF64vec(len);
      const PS::F64 sign = flap ? 1.0 : -1.0;
      base[axis] = 0.5 * len[axis] + ((Parameter::all_unit - 1) * Parameter::bond_leng + eps) * sign;
      nv[axis] = -sign;
      SetAmphilPartPos(base, nv, idx);

      if (asym) {
        const PS::U32 amp_id = idx / Parameter::all_unit;
        if (((amp_id % 100) == rat) || ((amp_id % 100) == 0 && (amp_id != 0)))
          flap ^= true;
      } else {
        flap ^= true;
      }
    }
  }

  bool MakeSphLine(const PS::F64 the,
                   const PS::F64vec& center,
                   const PS::F64 rad,
                   const PS::F64 sign,
                   const PS::U32 end_idx,
                   PS::U32& idx) {
    const PS::F64 q_thick = 0.5 * Parameter::all_unit * Parameter::bond_leng;
    const PS::F64 prj_rad = rad * std::sin(the);

    PS::F64 d_phi = lip_len / prj_rad;
    const PS::S32 phi_elem = static_cast<PS::S32>(2.0 * M_PI / d_phi);
    d_phi = 2.0 * M_PI / phi_elem;

    for (PS::S32 i = 0; i < phi_elem; i++) {
      const PS::F64 phi = d_phi * i;
      const PS::F64vec nv(sign * std::sin(the) * std::cos(phi),
                          sign * std::sin(the) * std::sin(phi),
                          sign * std::cos(the));
      const PS::F64vec cent2arc = sign * nv * (rad - sign * q_thick) + center;
      SetAmphilPartPos(cent2arc, nv, idx);
      if (idx >= end_idx) return false;
    }
    return true;
  }

  void MakeLineForEachTheta(const PS::S32 elem,
                            const PS::F64 d_the,
                            const PS::F64 offset,
                            const PS::F64 rad,
                            const PS::F64 sign,
                            const PS::U32 end_idx,
                            PS::U32& idx) {
    if (idx >= end_idx) return;
    const PS::F64vec cent = 0.5 * box_leng;
    for(PS::S32 i = 0; i < elem; i++) {
      const PS::F64 the = d_the * (i + offset);
      const bool flag = MakeSphLine(the, cent, rad, sign, end_idx, idx);
      if (!flag) return;
    }
  }

  void MakeSphSheet(const PS::U32 num,
                    const PS::U32 offset) {
    const PS::F64 min_box_leng = box_leng.getMin();
    if (!(sph_rad >= 0.0 && sph_rad <= 0.5 * min_box_leng)) {
      std::cerr << "sphere radius is too large.\n";
      std::cerr << "sph_rad: " << sph_rad << "min_box_leng " << min_box_leng << std::endl;
      std::exit(1);
    }

    const PS::F64 up_the = upper_the * M_PI / 180.0; //NOTE: the unit of upper_theta is not radian.

    const PS::F64 q_thick = 0.5 * Parameter::all_unit * Parameter::bond_leng;
    const PS::F64 out_rad = sph_rad + q_thick, in_rad  = sph_rad - q_thick;
    PS::F64 d_the_out = lip_len / out_rad, d_the_in = lip_len / in_rad;
    const PS::S32 out_the_elem = static_cast<PS::S32>(up_the / d_the_out), in_the_elem = static_cast<PS::S32>(up_the / d_the_in);
    d_the_out = up_the / out_the_elem;
    d_the_in  = up_the / in_the_elem;

    PS::U32 idx = offset;
    const PS::U32 end_idx = offset + num;
    MakeLineForEachTheta(out_the_elem, d_the_out, 0.0, out_rad, -1.0, end_idx, idx); //out
    MakeLineForEachTheta(in_the_elem, d_the_in, 0.0, in_rad, 1.0, end_idx, idx); //in

    const PS::U32 residue = end_idx - idx;
    if (residue > 0) {
      MakeLineForEachTheta(in_the_elem, d_the_in, 0.5, in_rad, 1.0, end_idx, idx);
      if (end_idx != idx) {
        std::cerr << "Too many particles are specified.\n";
        std::cerr << __FILE__ << " " << __LINE__ << std::endl;
        std::exit(1);
      }
    }
  }

  void MakeLineForEachAxis(const PS::S32 elem,
                           const PS::F64 rad,
                           const PS::F64 sign,
                           const PS::S32 axis,
                           const PS::U32 end_idx,
                           PS::U32& idx) {
    const PS::F64 h_pi = std::acos(0.0);

    PS::F64vec cent = 0.5 * box_leng;
    for (PS::S32 i = 0; i < elem; i++) {
      cent[axis] = i * lip_len;
      const bool flag = MakeSphLine(h_pi, cent, rad, sign, end_idx, idx);
      if (!flag) return;
    }
  }

  void MakeCylindSheet(const PS::S32 axis,
                       const PS::U32 num,
                       const PS::U32 offset) {
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

    const PS::S32 z_elem = static_cast<PS::S32>(box_leng[axis] / lip_len);
    PS::U32 idx = offset;
    const PS::U32 end_idx = offset + num;
    MakeLineForEachAxis(z_elem, out_rad, -1.0, axis, end_idx, idx); //out
    MakeLineForEachAxis(z_elem, in_rad, 1.0, axis, end_idx, idx);  //in

    const PS::U32 residue = end_idx - idx;
    if (residue > 0) {
      MakeLineForEachAxis(in_the_elem, out_rad, -1.0, axis, end_idx, idx);
      if (end_idx != idx) {
        std::cerr << "Too many particles are specified.\n";
        std::cerr << __FILE__ << " " << __LINE__ << std::endl;
        std::exit(1);
      }
    }
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

  void MakeCuboidSolvent(const PS::F64vec& up_pos,
                         const PS::F64vec& dw_pos,
                         const PS::U32 num,
                         const PS::U32 offset) {
    const PS::F64vec leng = up_pos - dw_pos;
    const PS::U32 beg_id = offset;
    const PS::U32 end_id = offset + num;
    for (PS::U32 id = beg_id; id < end_id; id++) {
      const PS::F64vec pos = GenrandF64vec(leng) + dw_pos;
      SetSolventPos(pos, id);
    }
  }

  void MakeSphericalSolvent(const bool is_in,
                            const PS::U32 num,
                            const PS::U32 offset,
                            const PS::U32 dim,
                            const PS::U32 axis) {
    const PS::U32 beg_id = offset;
    const PS::U32 end_id = offset + num;
    const PS::F64 h_thick = 0.5 * Parameter::all_unit * Parameter::bond_leng;

    PS::U32 idx = beg_id;
    while (idx < end_id) {
      const PS::F64vec pos = GenrandF64vec(box_leng);
      const PS::F64vec cent2pos = pos - 0.5 * box_leng;
      PS::F64 dist = cent2pos * cent2pos;
      if (dim == 2) {
        dist -= cent2pos[axis] * cent2pos[axis];
      }
      dist = std::sqrt(dist);
      if (is_in) {
        if (dist < sph_rad - h_thick) SetSolventPos(pos, idx++);
      } else {
        if (dist > sph_rad + h_thick) SetSolventPos(pos, idx++);
      }
    }
  }

#define MODE_EQ(str) strcasecmp(mode.c_str(), str) == 0
  void GenConfig() {
    if (MODE_EQ("flat")) {
      const PS::S32 axis = 1;
      MakeFlatSheet(box_leng, axis, false, amp_ptcl_num, 0);
      if (with_solvent) {
        PS::F64vec up_pos = box_leng, dw_pos(0.0, 0.0, 0.0);
        up_pos[axis] *= 0.5;
        MakeCuboidSolvent(up_pos, dw_pos, sol_num, amp_ptcl_num);
      }
    } else if (MODE_EQ("sphere")) {
      MakeSphSheet(amp_ptcl_num, 0);
      if (with_solvent) {
        MakeSphericalSolvent(true, sol_num, amp_ptcl_num, 3, 0);
      }
    } else if (MODE_EQ("cylind")) {
      const PS::S32 axis = 2; // Z axis
      MakeCylindSheet(axis, amp_ptcl_num, 0);
      if (with_solvent) {
        MakeSphericalSolvent(true, sol_num, amp_ptcl_num, 2, axis);
      }
    } else if (MODE_EQ("random")) {
      MakeRandomConfig(box_leng);
      if (with_solvent) {
        std::cerr << "with_solvent is not supported in this mode.\n";
        std::cerr << __FILE__ << " " << __LINE__ << std::endl;
        std::exit(1);
      }
    } else if (MODE_EQ("asym_flat")) {
      MakeFlatSheet(box_leng, 1, true, amp_ptcl_num, 0);
      if (with_solvent) {
        std::cerr << "with_solvent is not supported in this mode.\n";
        std::cerr << __FILE__ << " " << __LINE__ << std::endl;
        std::exit(1);
      }
    } else {
      if (!mode.empty())
        std::cerr << mode << ": Unknown mode\n";
      else
        std::cerr << "Mode is not specified\n.";
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;
      std::exit(1);
    }
  }

  void MatchingTagValues(std::map<std::string, std::vector<std::string>>& tag_val) {
    Parameter::Matching(&Tempera, std::string("Temperature"), tag_val, 1);
    Parameter::Matching(&(box_leng[0]), std::string("box_leng"), tag_val, 3);
    Parameter::Matching(&amp_num, "amp_num", tag_val, 1);
    amp_ptcl_num = Parameter::all_unit * amp_num;
    Parameter::Matching(&mode, "mode", tag_val, 1);

    if (tag_val.find("with_solvent") != tag_val.end()) {
      with_solvent = true;
    }

    if (MODE_EQ("sphere")) {
      Parameter::Matching(&sph_rad, "sph_rad", tag_val, 1);
      Parameter::Matching(&upper_the, "upper_the", tag_val, 1);
    }
    if (MODE_EQ("cylind")) {
      Parameter::Matching(&cyl_l, "cyl_l", tag_val, 1);
      Parameter::Matching(&cyl_r, "cyl_r", tag_val, 1);
    }
    if (MODE_EQ("asym_flat")) {
      Parameter::Matching(&in_out_rat, "in_out_rat", tag_val, 1);
    }
    if (with_solvent) {
      Parameter::Matching(&sol_num, "sol_num", tag_val, 1);
    }

    Parameter::Matching(&lip_len, "lip_len", tag_val, 1);
  }
#undef MODE_EQ

  void InitializeParticles() {
    for (auto& prtcl : prtcls) {
      prtcl.id = prtcl.prop = prtcl.amp_id = prtcl.unit = 0xffffffff;
      prtcl.pos.x = prtcl.pos.y = prtcl.pos.z = std::numeric_limits<PS::F64>::signaling_NaN();
      prtcl.vel.x = prtcl.vel.y = prtcl.vel.z = std::numeric_limits<PS::F64>::signaling_NaN();
      prtcl.vel_buf.x = prtcl.vel_buf.y = prtcl.vel_buf.z = std::numeric_limits<PS::F64>::signaling_NaN();
      prtcl.acc.x = prtcl.acc.y = prtcl.acc.z = std::numeric_limits<PS::F64>::signaling_NaN();
    }
  }

  void CheckParticleConfigIsValid() const {
    PS::F64 kin_temp = 0.0;
    for (const auto& prtcl : prtcls) {
      kin_temp += prtcl.vel * prtcl.vel;
      for (int j = 0; j < 3; j++) {
        if (!(prtcl.pos[j] <= box_leng[j] && prtcl.pos[j] >= 0.0)) {
          std::cerr << "There is a particle in outside range.\n";
          std::cerr << __FILE__ << " " << __LINE__ << std::endl;
          std::exit(1);
        }

        assert(std::isfinite(prtcl.vel_buf.x));
        assert(std::isfinite(prtcl.vel_buf.x));
        assert(std::isfinite(prtcl.vel_buf.x));

        assert(std::isfinite(prtcl.acc.x));
        assert(std::isfinite(prtcl.acc.x));
        assert(std::isfinite(prtcl.acc.x));
      }
    }

    kin_temp /= 3.0 * prtcls.size();
    if (fabs(kin_temp - Tempera) >= 1.0e-2) {
      std::cerr << "Please check kinetic temperature.\n";
      std::cerr << "Kinetic temperature is " << kin_temp << "K_BT.\n";
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;
      std::exit(1);
    }

    for (const auto& prtcl : prtcls) {
      if (!(prtcl.prop < Parameter::prop_num)) {
        std::cerr << "particle prop is invalid.\n";
        std::cerr << "prop is " << prtcl.prop << std::endl;
        std::exit(1);
      }
      if (prtcl.prop != Parameter::Solvent) {
        if (!(prtcl.amp_id < amp_num)) {
          std::cerr << "amp_id is invalid.\n";
          std::cerr << "amp_id is " << prtcl.amp_id << std::endl;
          std::exit(1);
        }
        if (!(prtcl.unit < Parameter::all_unit)) {
          std::cerr << "amphiphile unit is invalid\n";
          std::cerr << "unit is " << prtcl.unit << std::endl;
          std::exit(1);
        }
      } else {
        if (prtcl.amp_id != 0xffffffff) {
          std::cerr << "Solvent amp_id is not initialized.\n";
          std::exit(1);
        }
        if (prtcl.unit != 0xffffffff) {
          std::cerr << "Solvent unit is not initialized.\n";
          std::exit(1);
        }
      }
    }
  }

public:
  explicit ConfigMaker(const std::string& cdir_) {
    cdir = cdir_;
    box_leng.x = box_leng.y = box_leng.z = std::numeric_limits<PS::F64>::signaling_NaN();
  }
  ~ConfigMaker() {}

  void Initialize() {
    PS::MT::init_genrand(static_cast<unsigned long>(time(NULL)));
  }

  void LoadParam() {
    std::string fname = cdir + "/config_param.txt";
    std::ifstream fin(fname.c_str());
    if (fin) {
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
    PS::U32 all_n = amp_ptcl_num;
    if (with_solvent) all_n += sol_num;
    prtcls.resize(all_n);
    InitializeParticles();
    GenConfig();
    GenVeloc();
    CheckParticleConfigIsValid();
  }

  void DumpParticleConfig() {
    const std::string fname = cdir + "/init_config.xyz";
    FILE* fout = io_util::xfopen(fname.c_str(), "w");
    io_util::WriteXYZForm(&prtcls[0], prtcls.size(), 0, fout);
    fclose(fout);
  }
};

int main(int argc, char* argv[]) {
  if (argc != 2) {
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
