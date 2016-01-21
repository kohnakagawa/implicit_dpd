#pragma once

/*
  1. Random increase case
  2. 
 */

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
/*
  NOTE:
  1. Calculate the increase number of amphiphile.
  2. Do prefix-sum for all process.
  3. Calculate ids for each amphiphile molecule.
 */
class ChemManager {
public:
  
};

#else

//NOTE: PRNG is not thread safe.
template<class FP>
class ChemManager {
  static inline void Normalize(PS::F64vec& vec) {
    vec / std::sqrt(vec * vec);
  }

  static inline PS::F64 GenThermalVeloc() {
    return std::sqrt(-2.0 * Parameter::Tempera * std::log(PS::MT::genrand_real3() ) ) * std::cos(2.0 * M_PI * PS::MT::genrand_real3() );
  }

  static inline void ApplyPBC(PS::F64vec& pos,
			      const PS::F64vec& box_leng,
			      const PS::F64vec& ibox_leng) {
    pos.x -= std::floor(pos.x * ibox_leng.x) * box_leng.x;
    pos.y -= std::floor(pos.y * ibox_leng.y) * box_leng.y;
    pos.z -= std::floor(pos.z * ibox_leng.z) * box_leng.z;
  }

  template<PS::U32 part_num, PS::U32 unit_beg>
  void CreateAmpPart(PS::ParticleSystem<FP>& sys,
		     const PS::U32 prop, 
		     PS::U32 tpl_beg_id,
		     PS::U32& new_ptcl_id, const PS::U32 new_amp_id,
		     const PS::F64vec& tang_vec,
		     const PS::F64vec& box_leng, const PS::F64vec& ibox_leng,
		     PS::ReallocatableArray<PS::U32>& glob_topol) {
    FP new_ptcl;
    new_ptcl.prop = prop;
    
    for(PS::U32 i = 0; i < part_num; i++) {
      const PS::U32 base_pid = glob_topol[tpl_beg_id++];
      
      new_ptcl.id = new_ptcl_id;
      new_ptcl.amp_id = new_amp_id;
      new_ptcl.unit = unit_beg + i;
      new_ptcl.pos = sys[base_pid].pos + tang_vec;
      ApplyPBC(new_ptcl.pos, box_leng, ibox_leng);
      new_ptcl.delta_sumr = 0.0;
      
      new_ptcl.vel.x = GenThermalVeloc();
      new_ptcl.vel.y = GenThermalVeloc();
      new_ptcl.vel.z = GenThermalVeloc();

      new_ptcl.vel_buf = 0.0;
      new_ptcl.acc = 0.0;
      new_ptcl.press = 0.0;

      sys.AddNewParticle(new_ptcl);
      
      //Append topology list.
      glob_topol.pushBackNoCheck(new_ptcl_id);
      
      new_ptcl_id++;
    }
  }
  
public:
  ChemManager(const PS::U32 seed) {
    PS::MT::init_genrand(seed);
  }
  
  bool RandomChemEvent(PS::ParticleSystem<FP>& sys,
		       PS::ReallocatableArray<PS::U32>& glob_topol,
		       Parameter& param) {
    PS::U32 new_ptcl_id = sys.getNumberOfParticleLocal();
    PS::U32 new_amp_id = param.amp_num;
    const PS::S32 old_amp_num = param.amp_num;
    
    for(PS::S32 i = 0; i < old_amp_num; i++) {
      const PS::F64 rnd = PS::MT::genrand_real1();
      if(rnd <= param.p_thresld) {
	const PS::U32 head_id = glob_topol[Parameter::all_unit * i          ];
	const PS::U32 tail_id = glob_topol[Parameter::all_unit * (i + 1) - 1];
	PS::F64vec h2e = sys[tail_id].pos - sys[head_id].pos;
	ForceBonded<PS::ParticleSystem<FP> >::MinImage(h2e);
	PS::F64vec tang_vec(PS::MT::genrand_real1(), PS::MT::genrand_real1(), 0.0);
	tang_vec.z = -(h2e.x * tang_vec.x + h2e.y * tang_vec.y) / h2e.z;
	Normalize(tang_vec);
	tang_vec *= param.eps;

	const PS::U32 tpl_beg_id = Parameter::all_unit * i;
	CreateAmpPart<Parameter::head_unit, 0                   >(sys, Parameter::Hyphil, tpl_beg_id,
								  new_ptcl_id, new_amp_id, tang_vec,
								  param.box_leng, param.ibox_leng,
								  glob_topol);
	CreateAmpPart<Parameter::tail_unit, Parameter::head_unit>(sys, Parameter::Hyphob, tpl_beg_id + Parameter::head_unit,
								  new_ptcl_id, new_amp_id, tang_vec,
								  param.box_leng, param.ibox_leng,
								  glob_topol);
	param.amp_num++;
	new_amp_id++;
	
#ifdef DEBUG
	assert(param.amp_num == new_amp_id);
	assert(sys.getNumberOfParticleLocal() == new_ptcl_id);
#endif

	if(param.amp_num >= param.max_amp_num) {
	  std::cout << "Now param.amp_num >= param.max_amp_num.\n";
	  std::cout << "So, simulation will stop.\n";
	  return false;
	}
      }
    }
    return true;
  }

  bool CatalChemEvent(PS::ParticleSystem<FP>& sys,
		      Parameter& param) {
    //
  }
};

#endif
