#pragma once

#include <numeric>

template<class Tpsys>
struct ForceBonded {
  PS::ReallocatableArray<PS::U32> glob_topol;
  
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
  std::vector<PS::F64vec> buf_vir;
#endif

  ForceBonded(Tpsys& sys, const PS::U32 buf_size) {
    //NOTE: Bonded list construction is needed once when using OpenMP version.
    glob_topol.resizeNoInitialize(buf_size);
#ifdef DEBUG
    for(PS::U32 i = 0; i < buf_size; i++)
      glob_topol[i] = 0xffffffff;
#endif
    MakeGlobalBondedList(sys);

#ifdef DEBUG
    PS::U32 n = sys.getNumberOfParticleLocal();
    for(PS::U32 i = 0; i < buf_size; i++)
      assert(glob_topol[i] >= 0 && glob_topol[i] < n);
#endif

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
    buf_vir.resize(PS::Comm::getNumberOfThread(), 0.0);
#endif
    
  }
  ~ForceBonded() {}

  static inline void MinImage(PS::F64vec& drij) {
    drij.x -= Parameter::box_leng.x * std::round(drij.x * Parameter::ibox_leng.x);
    drij.y -= Parameter::box_leng.y * std::round(drij.y * Parameter::ibox_leng.y);
    drij.z -= Parameter::box_leng.z * std::round(drij.z * Parameter::ibox_leng.z);
  }

  void MakeGlobalBondedList(const Tpsys& sys) {
    PS::U32 n = sys.getNumberOfParticleLocal();
    
    for(PS::U32 i = 0; i < n; i++) {
      const PS::U32 aid = sys[i].amp_id;
      const PS::U32 unit = sys[i].unit;
      glob_topol[aid * Parameter::all_unit + unit] = i;
    }
  }

  void ExpandTopolBuffer(const PS::S32 n) {
    glob_topol.reserve(n);
  }

  //for intra cell
  static void StoreBondForceWithARLaw(const PS::F64vec&	__restrict dr,
				      const PS::F64&	__restrict inv_dr,
				      PS::F64vec&	__restrict d_vir,
				      PS::F64&		__restrict d_lap,
				      PS::F64vec*	__restrict F)
  {
    const PS::F64 cf_bond = Parameter::cf_spring<Parameter::bond_leng != 0.0>(inv_dr);
    
    const PS::F64vec Fbond(cf_bond * dr.x, cf_bond * dr.y, cf_bond * dr.z);

    //NOTE: The value of virial is twice.
    d_vir.x += 2.0 * dr.x * Fbond.x;
    d_vir.y += 2.0 * dr.y * Fbond.y;
    d_vir.z += 2.0 * dr.z * Fbond.z;

    //NOTE: The value of lap is twice.
    d_lap += 2.0 * Parameter::cf_s * (6.0 * Parameter::ibond - 4.0 * inv_dr);
    
    F[0] -= Fbond;
    F[1] += Fbond;
  }
  
  static void StoreBendForceWithARLaw(const PS::F64vec*	__restrict dr,
				      const PS::F64*	__restrict inv_dr,
				      const PS::F64*	__restrict dist,
				      PS::F64vec&	__restrict d_vir,
				      PS::F64&		__restrict d_lap,
				      PS::F64vec*	__restrict F)
  {
    const PS::F64	inv_dr_prod	= inv_dr[0] * inv_dr[1];
    const PS::F64	inv_dist[2]	= { inv_dr[0] * inv_dr[0],
					    inv_dr[1] * inv_dr[1] };
    const PS::F64	in_prod		= dr[0] * dr[1];
    const PS::F64	cf_bd		= Parameter::cf_b * inv_dr_prod;
    const PS::F64       cf_crs[2]	= { in_prod * inv_dist[0],
					    in_prod * inv_dist[1] };
    
    const PS::F64vec Ftb0(cf_bd * (dr[1].x - cf_crs[0] * dr[0].x),
			  cf_bd * (dr[1].y - cf_crs[0] * dr[0].y),
			  cf_bd * (dr[1].z - cf_crs[0] * dr[0].z));
    const PS::F64vec Ftb1(cf_bd * (dr[0].x - cf_crs[1] * dr[1].x),
			  cf_bd * (dr[0].y - cf_crs[1] * dr[1].y),
			  cf_bd * (dr[0].z - cf_crs[1] * dr[1].z));
    
    //NOTE: The value of virial is twice.
    d_vir.x += 2.0 * (dr[0].x * Ftb0.x + dr[1].x * Ftb1.x);
    d_vir.y += 2.0 * (dr[0].y * Ftb0.y + dr[1].y * Ftb1.y);
    d_vir.z += 2.0 * (dr[0].z * Ftb0.z + dr[1].z * Ftb1.z);

    //NOTE: The value of lap is twice.
    d_lap += 2.0 * 2.0 * cf_bd * inv_dist[0] * inv_dist[1] * ( in_prod * ( in_prod + 2.0 * (dist[0] + dist[1]) ) + dist[0] * dist[1]);
    
    F[0] -= Ftb0;
    F[1] += Ftb0 - Ftb1;
    F[2] += Ftb1;
  }

  //ASSUME: bond_n >= 3.
  template<PS::U32 bond_n>
  void CalcBondBendGlobalCell(Tpsys&		__restrict sys,
			      PS::U32		           beg_bond_id,
			      PS::F64vec&	__restrict d_vir,
			      PS::F64&		__restrict d_lap)
  {
    PS::F64vec Fbb[bond_n], pos_buf[bond_n], dr[bond_n - 1];
    PS::F64  dist2[bond_n - 1], inv_dr[bond_n - 1];

    pos_buf[0] = sys[ glob_topol[beg_bond_id    ] ].pos;
    pos_buf[1] = sys[ glob_topol[beg_bond_id + 1] ].pos;
    
    dr[0] = pos_buf[1] - pos_buf[0];
    MinImage(dr[0]);
    dist2[0] = dr[0] * dr[0];
    inv_dr[0] = 1.0 / std::sqrt(dist2[0]);
    
    StoreBondForceWithARLaw(dr[0], inv_dr[0], d_vir, d_lap, &Fbb[0]);

#pragma unroll
    for(PS::U32 unit = 2; unit < bond_n; unit++) {
      pos_buf[unit] = sys[ glob_topol[beg_bond_id + unit] ].pos;
      dr[unit - 1] = pos_buf[unit] - pos_buf[unit - 1];
      MinImage(dr[unit - 1]);
      dist2[unit - 1] = dr[unit - 1] * dr[unit - 1];
      inv_dr[unit - 1] = 1.0 / std::sqrt(dist2[unit - 1]);
      
      StoreBondForceWithARLaw(dr[unit - 1], inv_dr[unit - 1], d_vir, d_lap, &Fbb[unit - 1]);
      StoreBendForceWithARLaw(&dr[unit - 2], &inv_dr[unit - 2], dist2, d_vir, d_lap, &Fbb[unit - 2]);
    }

    //Store the sum of force.
#pragma unroll
    for(PS::U32 unit = 0; unit < bond_n; unit++)
      sys[ glob_topol[beg_bond_id + unit] ].acc += Fbb[unit];
  }

  void CalcListedForce(Tpsys& sys, PS::F64vec& bonded_vir) {
    //only intra cell
    const PS::U32 topol_num = glob_topol.size();
    
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
    
#pragma omp parallel
    {
      const PS::S32 tid = PS::Comm::getThreadNum();
      PS::F64vec d_vir(0.0); PS::F64 d_lap = 0.0;
#pragma omp for nowait
      for(PS::U32 i = 0; i < topol_num; i += Parameter::all_unit)
	CalcBondBendGlobalCell<Parameter::all_unit>(sys, i, d_vir, d_lap);
      buf_vir[tid] = d_vir;
    }
    bonded_vir = std::accumulate(buf_vir.cbegin(), buf_vir.cend(), PS::F64vec(0.0, 0.0, 0.0));
    
#else //no omp version

    PS::F64vec d_vir(0.0); PS::F64 d_lap = 0.0;
    for(PS::U32 i = 0; i < topol_num; i += Parameter::all_unit)
      CalcBondBendGlobalCell<Parameter::all_unit>(sys, i, d_vir, d_lap);
    bonded_vir = d_vir;
    
#endif
  }
  
  //no copy of virial
  void CalcListedForce(Tpsys& sys) {
    PS::F64vec buf(0.0, 0.0, 0.0);
    CalcListedForce(sys, buf);
  }
  
  //check topology
  void CheckBondedTopology(const Tpsys& sys) {
    PS::U32 n = sys.getNumberOfParticleLocal();
    for(PS::U32 i = 0; i < n; i++) {
      const PS::U32 aid = sys[i].amp_id;
      const PS::U32 unit = sys[i].unit;
      assert(glob_topol[aid * Parameter::all_unit + unit] == i);
    }
  }
  
  //for debug
  void CalcListedForceWithCheck(Tpsys& sys, PS::F64vec& bonded_vir) {
    CheckBondedTopology(sys);
    CalcListedForce(sys, bonded_vir);
  }
};


template<class Tpsys>
struct ForceBondedMPI {
  struct ampid2idx {
    PS::U32 amp_id, unit, idx;
    inline PS::U32 getKey() { return Parameter::all_unit * amp_id + unit; }
  };

  PS::U32 cmplt_amp = 0, imcmplt_amp = 0;
  PS::ReallocatableArray<ampid2idx> ampid, ampid_buf;
  PS::ReallocatableArray<PS::U32> loc_topol_cmpl, loc_topol_imcmpl;
  PS::RadixSort<PS::U32> rsorter;

  explicit ForceBondedMPI(const int est_loc_amp) {
    ampid.resizeNoInitialize(est_loc_amp);
    ampid_buf.resizeNoInitialize(est_loc_amp);
    
    loc_topol_cmpl.resizeNoInitialize(est_loc_amp * Parameter::all_unit);
    loc_topol_imcmpl.resizeNoInitialize(est_loc_amp * Parameter::all_unit);
  }
  ~ForceBondedMPI() {
    
  }

  //for inter cell
  static void StoreBondForceNoARLaw(const PS::F64vec&	__restrict dr,
				    const PS::F64&	__restrict inv_dr,
				    PS::F64vec&		__restrict d_vir,
				    PS::F64&		__restrict d_lap,
				    PS::F64vec*		__restrict F,
				    const bool*		__restrict mask)
  {
    const double cf_bond = Parameter::cf_spring<Parameter::bond_leng != 0.0>(inv_dr);
    
    const PS::F64vec Fbond(cf_bond * dr.x, cf_bond * dr.y, cf_bond * dr.z);
    
    d_vir.x += dr.x * Fbond.x;
    d_vir.y += dr.y * Fbond.y;
    d_vir.z += dr.z * Fbond.z;

    d_lap += Parameter::cf_s * (6.0 * Parameter::ibond - 4.0 * inv_dr);
    
    if(mask[0]) F[0] -= Fbond;
    if(mask[1]) F[1] += Fbond;
  }

  static void StoreBendForceNoARLaw(const PS::F64vec*	__restrict dr,
				    const PS::F64*	__restrict inv_dr,
				    const PS::F64*	__restrict dist,
				    PS::F64vec&		__restrict d_vir,
				    PS::F64&		__restrict d_lap,
				    PS::F64vec*		__restrict F,
				    const bool*		__restrict mask)
  {
    const PS::F64	inv_dr_prod	= inv_dr[0] * inv_dr[1];
    const PS::F64	inv_dist[2]	= { inv_dr[0] * inv_dr[0],
					    inv_dr[1] * inv_dr[1] };
    const PS::F64	in_prod		= dr[0] * dr[1];
    const PS::F64	cf_bd		= Parameter::cf_b * inv_dr_prod;
    const PS::F64       cf_crs[2]	= { in_prod * inv_dist[0],
					    in_prod * inv_dist[1] };
    
    const PS::F64vec Ftb0(cf_bd * (dr[1].x - cf_crs[0] * dr[0].x),
			  cf_bd * (dr[1].y - cf_crs[0] * dr[0].y),
			  cf_bd * (dr[1].z - cf_crs[0] * dr[0].z));
    const PS::F64vec Ftb1(cf_bd * (dr[0].x - cf_crs[1] * dr[1].x),
			  cf_bd * (dr[0].y - cf_crs[1] * dr[1].y),
			  cf_bd * (dr[0].z - cf_crs[1] * dr[1].z));
    
    d_vir.x += dr[0].x * Ftb0.x + dr[1].x * Ftb1.x;
    d_vir.y += dr[0].y * Ftb0.y + dr[1].y * Ftb1.y;
    d_vir.z += dr[0].z * Ftb0.z + dr[1].z * Ftb1.z;

    d_lap += 2.0 * cf_bd * inv_dist[0] * inv_dist[1] * ( in_prod * ( in_prod + 2.0 * (dist[0] + dist[1]) ) + dist[0] * dist[1]);
    
    if(mask[0]) F[0] -= Ftb0;
    if(mask[1]) F[1] += Ftb0 - Ftb1;
    if(mask[2]) F[2] += Ftb1;    
  }

  //NOTE: minimum image convention is not needed.
  //ASSUME: bond_n >= 3.
  template<PS::U32 bond_n>
  void CalcBondBendLocalCell(Tpsys&		__restrict sys,
			     PS::U32		           beg_bond_id,
			     PS::F64vec&	__restrict d_vir,
			     PS::F64&		__restrict d_lap)
  {
    PS::F64vec Fbb[bond_n], pos_buf[bond_n], dr[bond_n - 1];
    PS::F64  dist2[bond_n - 1], inv_dr[bond_n - 1];

    const PS::U32 l_dst[] = { loc_topol_cmpl[beg_bond_id], loc_topol_cmpl[beg_bond_id + 1] };
    pos_buf[0] = sys[ l_dst[0] ].pos; pos_buf[1] = sys[ l_dst[1] ].pos; 
    
    dr[0] = pos_buf[1] - pos_buf[0];
    dist2[0] = dr[0] * dr[0];
    inv_dr[0] = 1.0 / std::sqrt(dist2[0]);
    
    ForceBonded<Tpsys>::StoreBondForceWithARLaw(dr[0], inv_dr[0], d_vir, d_lap, &Fbb[0]);

#pragma unroll
    for(PS::U32 unit = 2; unit < bond_n; unit++) {
      pos_buf[unit] = sys[ loc_topol_cmpl[beg_bond_id + unit] ].pos;
      dr[unit - 1] = pos_buf[unit] - pos_buf[unit - 1];
      dist2[unit - 1] = dr[unit - 1] * dr[unit - 1];
      inv_dr[unit - 1] = 1.0 / std::sqrt(dist2[unit - 1]);
      
      ForceBonded<Tpsys>::StoreBondForceWithARLaw(dr[unit - 1], inv_dr[unit - 1], d_vir, d_lap, &Fbb[unit - 1]);
      ForceBonded<Tpsys>::StoreBendForceWithARLaw(&dr[unit - 2], &inv_dr[unit - 2], dist2, d_vir, d_lap, &Fbb[unit - 2]);
    }

    //Store the sum of force.
#pragma unroll
    for(PS::U32 unit = 0; unit < bond_n; unit++)
      sys[ loc_topol_cmpl[beg_bond_id + unit] ].acc += Fbb[unit];
  }

  void MakeLocalBondedList(const Tpsys& sys) {
    PS::U32 n = sys.getNumberOfParticleLocal();
    
    for(PS::U32 i = 0; i < n; i++) {
      ampid[i].amp_id	= sys[i].amp_id;
      ampid[i].unit	= sys[i].unit;
      ampid[i].idx	= i;
    }
    rsorter.lsdSort(ampid, ampid_buf, 0, n);
    
    cmplt_amp = imcmplt_amp = 0;
    PS::U32 id_cmpl = 0, id_imcmpl = 0, cnt = 1, id_bef = ampid[0].amp_id, id_cur;
    for(PS::U32 i = 1; i < n; i++) {
      id_cur = ampid[i].amp_id;
      if(id_cur == id_bef) {
	cnt++;
      } else {
	if(cnt == Parameter::all_unit) {
	  for(PS::U32 j = 0; j < Parameter::all_unit; j++)
	    loc_topol_cmpl[id_cmpl++] = ampid[i - j - 1].idx;
	  cmplt_amp++;
	} else {
	  PS::U32 imcmpl_buf[Parameter::all_unit] = { 0xffffffff };
	  for(PS::U32 j = 0; j < cnt; j++)
	    imcmpl_buf[ ampid[i - j - 1].unit ] = ampid[i - j - 1].idx;
	  for(PS::U32 j = 0; j < Parameter::all_unit; j++)
	    loc_topol_imcmpl[id_imcmpl++] = imcmpl_buf[j];
	  imcmplt_amp++;
	}
	cnt = 1;
      }
      id_bef = id_cur;
    }
  }

  void MakeSurfaceBondedList(const Tpsys& sys) {
    
  }

  void CalcListedForce(Tpsys& sys) {
    //MPI
    //for intra cell
    MakeLocalBondedList(sys);
    PS::F64vec d_vir(0.0); PS::F64 d_lap = 0.0;
    for(PS::U32 i = 0; i < cmplt_amp; i++)
      CalcBondBendLocalCell<Parameter::all_unit>(sys, i * Parameter::all_unit, d_vir, d_lap);

    //for inter cell
    
  }
};

