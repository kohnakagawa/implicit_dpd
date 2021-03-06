#pragma once

template<class FP>
class ChemManager {
  PS::ReallocatableArray<PS::U32> target_id;
  PS::ReallocatableArray<PS::F64vec> core_poss_h_, h2t_vecs_;

  static inline void Normalize(PS::F64vec& vec) {
    vec /= std::sqrt(vec * vec);
  }

  static inline PS::F64vec CrossVecWithNormalize(const PS::F64vec& a0, const PS::F64vec a1) {
    PS::F64vec v = a0 ^ a1;
    Normalize(v);
    return v;
  }

  static inline PS::F64 GenThermalVeloc() {
    return std::sqrt(-2.0 * Parameter::Tempera * std::log(PS::MT::genrand_real3())) * std::cos(2.0 * M_PI * PS::MT::genrand_real3());
  }

  static inline void ApplyPBC(PS::F64vec& pos) {
    pos.x -= std::floor(pos.x * Parameter::ibox_leng.x) * Parameter::box_leng.x;
    pos.y -= std::floor(pos.y * Parameter::ibox_leng.y) * Parameter::box_leng.y;
    pos.z -= std::floor(pos.z * Parameter::ibox_leng.z) * Parameter::box_leng.z;
  }

  inline void AppendTopol(PS::ReallocatableArray<PS::U32>& glob_topol,
			  const PS::U32 part_num,
			  const PS::U32 new_ptcl_id) {
    for (PS::U32 i = 0; i < part_num; i++)
      glob_topol.pushBackNoCheck(new_ptcl_id + i);
  }

  template<PS::U32 part_num, PS::U32 unit_beg>
  void CreateAmpPart(PS::ParticleSystem<FP>& sys,
                     const PS::U32 prop,
                     PS::U32& new_ptcl_id,
                     const PS::U32 new_amp_id,
                     const PS::F64vec& copied_pos,
                     const PS::F64vec& h2e) {
    FP new_ptcl;
    new_ptcl.prop = prop;
    new_ptcl.amp_id = new_amp_id;

    for (PS::U32 i = 0; i < part_num; i++) {
      new_ptcl.id   = new_ptcl_id;
      new_ptcl.unit = unit_beg + i;
      new_ptcl.pos  = copied_pos + h2e * (unit_beg + i);
      ApplyPBC(new_ptcl.pos);
      new_ptcl.delta_sumr = 0.0;

      new_ptcl.vel.x = GenThermalVeloc();
      new_ptcl.vel.y = GenThermalVeloc();
      new_ptcl.vel.z = GenThermalVeloc();

      new_ptcl.vel_buf = 0.0;
      new_ptcl.acc = 0.0;
      new_ptcl.press = 0.0;

      new_ptcl.density.fill(0.0);

      sys.AddNewParticle(new_ptcl);

      new_ptcl_id++;
    }
  }

  void DetermineTargetId(const PS::U32 amp_num,
                         const PS::ReallocatableArray<PS::U32>& topol,
                         const PS::U32 loc_num,
                         Parameter& param) {
    for (PS::U32 i = 0; i < amp_num; i++) {
      const PS::F64 rnd = PS::MT::genrand_real1();
      const PS::U32 head_id = topol[Parameter::all_unit * i];
      if ((rnd <= param.p_thresld) && (head_id < loc_num)) { // only consider real particle
        const PS::U32 tail_id = topol[Parameter::all_unit * i + 2];
        target_id.push_back(head_id);
        target_id.push_back(tail_id);
      }
    }
  }

  void GetCorePos(PS::ParticleSystem<FP>& sys,
                  const PS::U32 loc_num,
                  const PS::U32 core_id,
                  PS::F64vec& core_pos_h,
                  PS::F64vec& core_pos_t,
                  PS::F64vec& core_pos_cm,
                  const Parameter& param) {
    core_pos_h = -1.0;
    core_pos_t = -1.0;
    core_pos_cm = -1.0;
    for (PS::U32 i = 0; i < loc_num; i++) {
      if (sys[i].id == param.core_ptcl_id[core_id]) {
        core_pos_h  = sys[i].pos;
        core_pos_t  = sys[i].nei_cm_pos[Parameter::Hyphob];
        core_pos_cm = sys[i].nei_cm_pos[Parameter::Hyphil];
        break;
      }
    }

    PS::S32 rank = PS::Comm::getRank();
    if ((core_pos_h.x < 0.0) &&
        (core_pos_h.y < 0.0) &&
        (core_pos_h.z < 0.0)) {
      rank = -1;
    }

    const PS::S32 root_rank = PS::Comm::getMaxValue(rank);
    PS::Comm::broadcast(&core_pos_h, 1, root_rank);
    PS::Comm::broadcast(&core_pos_t, 1, root_rank);
    PS::Comm::broadcast(&core_pos_cm, 1, root_rank);
  }

  void DetermineTargetIdWithCorePos(PS::ParticleSystem<FP>& sys,
                                    const PS::U32 amp_num,
                                    const PS::ReallocatableArray<PS::U32>& topol,
                                    const PS::U32 loc_num,
                                    Parameter& param) {
    PS::F64vec core_pos_h, core_pos_t, core_pos_cm;
    const PS::U32 num_core_amp_id = param.core_amp_id().size();

    for (PS::U32 i = 0; i < num_core_amp_id; i++) {
      GetCorePos(sys, loc_num, i, core_pos_h, core_pos_t, core_pos_cm, param);

      // first calculate membrane normal vector
      PS::F64vec h2t = core_pos_t - core_pos_h;
      ForceBonded<PS::ParticleSystem<FP> >::MinImage(h2t);
      Normalize(h2t);

      // then reduce protorusion height fluct.
      PS::F64vec core2cmpos = core_pos_cm - core_pos_h;
      ForceBonded<PS::ParticleSystem<FP> >::MinImage(core2cmpos);
      const PS::F64 prj_cf = core2cmpos * h2t;
      core_pos_h += prj_cf * h2t;
      ApplyPBC(core_pos_h);

      // calc target id
      for (PS::U32 j = 0; j < amp_num; j++) {
        const PS::U32 head_id = topol[Parameter::all_unit * j];
        if (head_id < loc_num) { // only consider real particle
          PS::F64 rnd = 2.0;
          PS::F64vec core2ptcl = sys[head_id].pos - core_pos_h;
          ForceBonded<PS::ParticleSystem<FP> >::MinImage(core2ptcl);

          const PS::F64 depth = h2t * core2ptcl;
          const PS::F64vec tang_vec = core2ptcl - depth * h2t;
          const PS::F64 rad = std::sqrt(tang_vec * tang_vec);

          const PS::F64 rad_thresld = param.influ_grd * depth + param.influ_rad;

          if ((rad <= rad_thresld) && (depth >= -param.influ_hei) && (depth <= param.influ_dep))
            rnd = PS::MT::genrand_real1();

          if (rnd <= param.p_thresld) {
            const PS::U32 tail_id = topol[Parameter::all_unit * j + 2];
            target_id.push_back(head_id);
            target_id.push_back(tail_id);
          }
        }
      }
    }
  }

  PS::F64vec GetMembraneTangRandomVec(const PS::F64vec& h2e) const {
    const PS::F64 theta = 2.0 * M_PI * PS::MT::genrand_real1();
    const PS::F64vec base(std::cos(theta), std::sin(theta), 0.0);
    const PS::F64 base_dot_h2e = base * h2e;
    PS::F64vec tang_vec = base - h2e * base_dot_h2e;
    Normalize(tang_vec);
    return tang_vec;
  }

public:
  explicit ChemManager(const PS::U32 seed) {
    // NOTE: This PRNG is not thread safe.
    PS::MT::init_genrand(seed);
    target_id.resizeNoInitialize(10000);
  }
  ~ChemManager() {}

  const PS::ReallocatableArray<PS::F64vec>& h2t_vecs() const {
    return h2t_vecs_;
  }

  const PS::ReallocatableArray<PS::F64vec>& core_poss_h() const {
    return core_poss_h_;
  }

  // MPI version
  template<class Pepj>
  bool RandomChemEvent(PS::ParticleSystem<FP>& sys,
                       const PS::ReallocatableArray<Pepj>& epj_org,
                       const PS::ReallocatableArray<PS::U32>& loc_topol_cmpl,
                       const PS::ReallocatableArray<PS::U32>& loc_topol_imcmpl,
                       const PS::U32 cmplt_amp,
                       const PS::U32 imcmplt_amp,
                       Parameter& param) {
    target_id.clearSize();

    const PS::U32 loc_num = sys.getNumberOfParticleLocal();
#ifdef LOCAL_CHEM_EVENT
    DetermineTargetIdWithCorePos(sys, cmplt_amp, loc_topol_cmpl, loc_num, param);
    DetermineTargetIdWithCorePos(sys, imcmplt_amp, loc_topol_imcmpl, loc_num, param);
#else
    DetermineTargetId(cmplt_amp, loc_topol_cmpl, loc_num, param);
    DetermineTargetId(imcmplt_amp, loc_topol_imcmpl, loc_num, param);
#endif

    const PS::U32 incr_num = target_id.size() / 2;
    PS::U32 offset = 0;      // amphiphile offset number
    PS::Comm::exScan(&incr_num, &offset, 1);
    PS::U32 new_ptcl_id = sys.getNumberOfParticleGlobal() + offset * Parameter::all_unit;
    offset += param.amp_num; // add global offset

    for (PS::U32 i = 0; i < incr_num; i++) {
      PS::F64vec h2e = epj_org[target_id[2 * i + 1]].pos - epj_org[target_id[2 * i]].pos;
      ForceBonded<PS::ParticleSystem<FP> >::MinImage(h2e);
      Normalize(h2e);
      const PS::F64vec tang_vec = GetMembraneTangRandomVec(h2e) * param.eps;

      h2e *= Parameter::bond_leng;
      const PS::F64vec copied_pos = tang_vec + epj_org[target_id[2 * i]].pos;
      const PS::U32 new_amp_id = offset + i;
      CreateAmpPart<Parameter::head_unit, 0                   >(sys, Parameter::Hyphil, new_ptcl_id,
                                                                new_amp_id, copied_pos, h2e);
      CreateAmpPart<Parameter::tail_unit, Parameter::head_unit>(sys, Parameter::Hyphob, new_ptcl_id,
                                                                new_amp_id, copied_pos, h2e);
    }

    param.amp_num += PS::Comm::getSum(incr_num);

    if (param.amp_num >= param.max_amp_num) {
      std::cout << "Now param.amp_num >= param.max_amp_num.\n";
      std::cout << "So, simulation will stop.\n";
      return false;
    } else {
      return true;
    }
  }

  // single core ver
  bool RandomChemEvent(PS::ParticleSystem<FP>& sys,
                       PS::ReallocatableArray<PS::U32>& glob_topol,
                       Parameter& param) {
    target_id.clearSize();

    const PS::U32 glb_num = sys.getNumberOfParticleGlobal();
#ifdef LOCAL_CHEM_EVENT
    DetermineTargetIdWithCorePos(sys, param.amp_num, glob_topol, glb_num, param);
#else
    DetermineTargetId(param.amp_num, glob_topol, glb_num, param);
#endif

    const PS::U32 incr_num = target_id.size() / 2;
    PS::U32 new_ptcl_id = glb_num;
    for (PS::U32 i = 0; i < incr_num; i++) {
      PS::F64vec h2e = sys[target_id[2 * i + 1]].pos - sys[target_id[2 * i]].pos;
      ForceBonded<PS::ParticleSystem<FP>>::MinImage(h2e);
      Normalize(h2e);
      const PS::F64vec tang_vec = GetMembraneTangRandomVec(h2e) * param.eps;

      h2e *= Parameter::bond_leng;
      const PS::F64vec copied_pos = tang_vec + sys[target_id[2 * i]].pos;
      const PS::U32 new_amp_id = param.amp_num + i;
      CreateAmpPart<Parameter::head_unit, 0                   >(sys, Parameter::Hyphil, new_ptcl_id,
                                                                new_amp_id, copied_pos, h2e);
      CreateAmpPart<Parameter::tail_unit, Parameter::head_unit>(sys, Parameter::Hyphob, new_ptcl_id,
                                                                new_amp_id, copied_pos, h2e);
    }
    param.amp_num += incr_num;

    if (param.amp_num >= param.max_amp_num) {
      std::cout << "Now param.amp_num >= param.max_amp_num.\n";
      std::cout << "So, simulation will stop.\n";
      return false;
    } else {
      return true;
    }
  }
};
