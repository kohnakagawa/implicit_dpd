#pragma once

template<class FP>
class ChemManager {
  PS::ReallocatableArray<PS::U32> target_id;

  static inline void Normalize(PS::F64vec& vec) {
    vec / std::sqrt(vec * vec);
  }

  static inline PS::F64 GenThermalVeloc() {
    return std::sqrt(-2.0 * Parameter::Tempera * std::log(PS::MT::genrand_real3() ) ) * std::cos(2.0 * M_PI * PS::MT::genrand_real3() );
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
      if ( (rnd <= param.p_thresld) && (head_id < loc_num) ) {
	const PS::U32 tail_id = topol[Parameter::all_unit * i + 2];
	target_id.push_back(head_id);
	target_id.push_back(tail_id);
      }
    }
  }
  
public:
  ChemManager(const PS::U32 seed) {
    //NOTE: This PRNG is not thread safe.
    PS::MT::init_genrand(seed);
    target_id.resizeNoInitialize(10000);
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
    DetermineTargetId(cmplt_amp, loc_topol_cmpl, loc_num, param);
    DetermineTargetId(imcmplt_amp, loc_topol_imcmpl, loc_num, param);
    const PS::U32 incr_num = target_id.size() / 2;
    PS::U32 offset = 0;      // amphiphile offset number
    PS::Comm::exScan(&incr_num, &offset, 1);
    PS::U32 new_ptcl_id = sys.getNumberOfParticleGlobal() + offset * Parameter::all_unit;
    offset += param.amp_num; // add global offset
    
    for(PS::U32 i = 0; i < incr_num; i++) {
      PS::F64vec h2e = epj_org[target_id[2 * i]].pos - epj_org[target_id[2 * i + 1]].pos;
      ForceBonded<PS::ParticleSystem<FP> >::MinImage(h2e);
      PS::F64vec tang_vec(PS::MT::genrand_real1(), PS::MT::genrand_real1(), 0.0);
      tang_vec.z = -(h2e.x * tang_vec.x + h2e.y * tang_vec.y) / h2e.z;
      Normalize(tang_vec);
      Normalize(h2e);
      tang_vec *= param.eps;
      h2e *= Parameter::bond_leng;

      const PS::F64vec copied_pos = tang_vec + epj_org[target_id[2 * i]].pos;
      const PS::U32 new_amp_id = offset + i;
      CreateAmpPart<Parameter::head_unit, 0                   >(sys, Parameter::Hyphil, new_ptcl_id,
								new_amp_id, copied_pos, h2e);
      CreateAmpPart<Parameter::tail_unit, Parameter::head_unit>(sys, Parameter::Hyphob, new_ptcl_id,
								new_amp_id, copied_pos, h2e);
    }

    param.amp_num += PS::Comm::getSum(incr_num);

    if(param.amp_num >= param.max_amp_num) {
      std::cout << "Now param.amp_num >= param.max_amp_num.\n";
      std::cout << "So, simulation will stop.\n";
      return false;
    } else {
      return true;
    }
  }
  
  // normal ver
  bool RandomChemEvent(PS::ParticleSystem<FP>& sys,
		       PS::ReallocatableArray<PS::U32>& glob_topol,
		       Parameter& param) {
    PS::U32 new_ptcl_id = sys.getNumberOfParticleLocal();
    PS::U32 new_amp_id = param.amp_num;
    const PS::U32 old_amp_num = param.amp_num;
    
    for (PS::U32 i = 0; i < old_amp_num; i++) {
      const PS::F64 rnd = PS::MT::genrand_real1();
      if (rnd <= param.p_thresld) {
	const PS::U32 head_id = glob_topol[Parameter::all_unit * i    ];
	const PS::U32 tail_id = glob_topol[Parameter::all_unit * i + 2];
	PS::F64vec h2e = sys[tail_id].pos - sys[head_id].pos;
	ForceBonded<PS::ParticleSystem<FP> >::MinImage(h2e);
	PS::F64vec tang_vec(PS::MT::genrand_real1(), PS::MT::genrand_real1(), 0.0);
	tang_vec.z = -(h2e.x * tang_vec.x + h2e.y * tang_vec.y) / h2e.z;
	Normalize(tang_vec);
	Normalize(h2e);
	tang_vec *= param.eps;
	h2e *= Parameter::bond_leng;
	
	const PS::F64vec copied_pos = tang_vec + sys[head_id].pos;
	AppendTopol(glob_topol, Parameter::head_unit, new_ptcl_id);
	CreateAmpPart<Parameter::head_unit, 0                   >(sys, Parameter::Hyphil, new_ptcl_id,
								  new_amp_id, copied_pos, h2e);
	AppendTopol(glob_topol, Parameter::tail_unit, new_ptcl_id);
	CreateAmpPart<Parameter::tail_unit, Parameter::head_unit>(sys, Parameter::Hyphob, new_ptcl_id,
								  new_amp_id, copied_pos, h2e);
	param.amp_num++;
	new_amp_id++;
	
#ifdef DEBUG
	assert(param.amp_num == new_amp_id);
	assert(sys.getNumberOfParticleLocal() == new_ptcl_id);
#endif

	if (param.amp_num >= param.max_amp_num) {
	  std::cout << "Now param.amp_num >= param.max_amp_num.\n";
	  std::cout << "So, simulation will stop.\n";
	  return false;
	}
      }
    }
    return true;
  }
  /*  bool CatalChemEvent(PS::ParticleSystem<FP>& sys,
		      Parameter& param) {
    //
    }*/
};
