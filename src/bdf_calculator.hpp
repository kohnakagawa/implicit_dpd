#pragma once

#include <numeric>
#include <algorithm>
#include <parallel/algorithm>

PS::F64 Parameter::cf_s;
PS::F64 Parameter::cf_b;

template<class Tpsys>
struct ForceBonded {
  PS::ReallocatableArray<PS::U32> glob_topol_;
  PS::ReallocatableArray<PS::F64> cf_bend;

  PS::ReallocatableArray<PS::U32>& glob_topol() {
    return glob_topol_;
  }

  // ASSUME: there are only amphiphile molecules in this system.
  void AppendCfBend() {
    const PS::S32 incr_num = glob_topol_.size() / Parameter::all_unit - cf_bend.size();
    for (PS::S32 i = 0; i < incr_num; i++) {
      cf_bend.push_back(Parameter::cf_b);
    }
  }
  
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
  std::vector<PS::F64vec> buf_vir;
#endif

  ForceBonded(Tpsys& sys,
	      const PS::U32 buf_size,
	      const PS::U32 amp_num,
	      const std::vector<PS::U32>& core_amp_id) {
    //NOTE: Bonded list construction is needed once when using OpenMP version.
    glob_topol_.resizeNoInitialize(buf_size);
    cf_bend.resizeNoInitialize(amp_num);
#ifdef DEBUG
    for (PS::U32 i = 0; i < buf_size; i++) glob_topol_[i] = 0xffffffff;
#endif
    MakeGlobalBondedList(sys, core_amp_id);

#ifdef DEBUG
    for (PS::U32 i = 0; i < buf_size; i++) assert(glob_topol_[i] < buf_size);
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

  void MakeGlobalBondedList(const Tpsys& sys, const std::vector<PS::U32>& core_amp_id) {
    PS::U32 n = sys.getNumberOfParticleLocal();
    PS::S32 cnt = 0;
    
    for (PS::U32 i = 0; i < n; i++) {
      const PS::U32 aid = sys[i].amp_id;
      const PS::U32 unit = sys[i].unit;
      if ((aid != 0xffffffff) && (unit != 0xffffffff)) {
	glob_topol_[aid * Parameter::all_unit + unit] = i;
	cnt++;
      } else {
	// This particle must be solvent.
	assert(sys[i].prop == Parameter::Solvent);
      }
    }
    
    const PS::U32 num_amp = cf_bend.size();
    const PS::U32 num_core_amp_id = core_amp_id.size();
    for (PS::U32 i = 0; i < num_amp; i++) {
      cf_bend[i] = Parameter::cf_b;
      for (PS::U32 j = 0; j < num_core_amp_id; j++) {
	if (i == core_amp_id[j]) {
	  cf_bend[i] = Parameter::cf_b_rigid;
	}
      }
    }
    
    CHECK_EQ(cnt, glob_topol_.size());
  }

  void ExpandTopolBuffer(const PS::U32 n) {
    glob_topol_.reserve(n);
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
				      PS::F64vec*	__restrict F,
				      const PS::F64     cf_b)
  {
    const PS::F64	inv_dr_prod	= inv_dr[0] * inv_dr[1];
    const PS::F64	inv_dist[2]	= { inv_dr[0] * inv_dr[0],
					    inv_dr[1] * inv_dr[1] };
    const PS::F64	in_prod		= dr[0] * dr[1];
    const PS::F64	cf_bd		= cf_b * inv_dr_prod;
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
    lap_conf += 2.0 * 2.0 * cf_bd * (in_prod * (2.0 * (inv_dist[0] + inv_dist[1]) + in_prod * inv_dr_prod * inv_dr_prod) + 1.0);
    
    F[0] -= Ftb0;
    F[1] += Ftb0 - Ftb1;
    F[2] += Ftb1;
  }

  //ASSUME: bond_n >= 3.
  template<PS::U32 bond_n>
  void CalcBondBendGlobalCell(Tpsys&		__restrict sys,
			      const PS::U32        beg_bond_id,
			      const PS::F64        cf_b,
			      PS::F64vec&	__restrict d_vir,
			      PS::F64&		__restrict d_lap)
  {
    PS::F64vec Fbb[bond_n], pos_buf[bond_n], dr[bond_n - 1];
    PS::F64  dist2[bond_n - 1], inv_dr[bond_n - 1];

    const PS::U32* l_dst = &(glob_topol_[beg_bond_id]);

    pos_buf[0] = sys[ l_dst[0] ].pos;
    pos_buf[1] = sys[ l_dst[1] ].pos;

    dr[0] = pos_buf[1] - pos_buf[0];
    MinImage(dr[0]);
    dist2[0] = dr[0] * dr[0];
    inv_dr[0] = 1.0 / std::sqrt(dist2[0]);
    
    StoreBondForceWithARLaw(dr[0], inv_dr[0], d_vir, d_lap, &Fbb[0]);

#pragma unroll
    for (PS::U32 unit = 2; unit < bond_n; unit++) {
      pos_buf[unit] = sys[ l_dst[unit] ].pos;
      dr[unit - 1] = pos_buf[unit] - pos_buf[unit - 1];
      MinImage(dr[unit - 1]);
      dist2[unit - 1] = dr[unit - 1] * dr[unit - 1];
      inv_dr[unit - 1] = 1.0 / std::sqrt(dist2[unit - 1]);
      
      StoreBondForceWithARLaw(dr[unit - 1], inv_dr[unit - 1], d_vir, d_lap, &Fbb[unit - 1]);
      StoreBendForceWithARLaw(&dr[unit - 2], &inv_dr[unit - 2], dist2, d_vir, d_lap, &Fbb[unit - 2], cf_b);
    }

    //Store the sum of force.
#pragma unroll
    for (PS::U32 unit = 0; unit < bond_n; unit++)
      sys[ glob_topol_[beg_bond_id + unit] ].acc += Fbb[unit];
  }

  void CalcListedForce(Tpsys& sys, PS::F64vec& bonded_vir) {
    // only intra cell
    AppendCfBend();

    const PS::U32 amp_num = cf_bend.size();
    
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
    
#pragma omp parallel
    {
      const PS::S32 tid = PS::Comm::getThreadNum();
      PS::F64vec d_vir(0.0); PS::F64 d_lap = 0.0;
#pragma omp for nowait
      for (PS::U32 i = 0; i < amp_num; i++) {
	CalcBondBendGlobalCell<Parameter::all_unit>(sys, Parameter::all_unit * i, cf_bend[i], d_vir, d_lap);
      }
      buf_vir[tid] = d_vir;
    }
    bonded_vir = std::accumulate(buf_vir.cbegin(), buf_vir.cend(), PS::F64vec(0.0, 0.0, 0.0));
    
#else //no omp version

    PS::F64vec d_vir(0.0); PS::F64 d_lap = 0.0;
    for (PS::U32 i = 0; i < amp_num; i++) {
      CalcBondBendGlobalCell<Parameter::all_unit>(sys, Parameter::all_unit * i, cf_bend[i], d_vir, d_lap);
    }
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
    for (PS::U32 i = 0; i < n; i++) {
      const PS::U32 aid = sys[i].amp_id;
      const PS::U32 unit = sys[i].unit;
      assert(glob_topol_[aid * Parameter::all_unit + unit] == i);
    }
  }
  
  //for debug
  void CalcListedForceWithCheck(Tpsys& sys, PS::F64vec& bonded_vir) {
    CheckBondedTopology(sys);
    CalcListedForce(sys, bonded_vir);
  }
};

template<class Tpsys, class Pepj>
struct ForceBondedMPI {
  struct ampid2idx {
    PS::U32 amp_id, unit, idx, key;
    bool is_real;
    inline void setKey() {
      if ((amp_id != 0xffffffff) && (unit != 0xffffffff))
    	key = Parameter::all_unit * amp_id + unit;
      else
    	key = 0xffffffff;
    }
  };

  PS::U32 cmplt_amp_ = 0, imcmplt_amp_ = 0;
  PS::ReallocatableArray<ampid2idx> ampid, ampid_buf;
  PS::ReallocatableArray<PS::U32> loc_topol_cmpl_, loc_topol_imcmpl_;
  PS::ReallocatableArray<bool> is_real_surf; //real particle or phantom on surface
  PS::ReallocatableArray<PS::F64> cf_bend_cmpl, cf_bend_imcmpl;
  
  PS::U32 cmplt_amp()   const { return cmplt_amp_; }
  PS::U32 imcmplt_amp() const { return imcmplt_amp_; }
  const PS::ReallocatableArray<PS::U32>& loc_topol_cmpl()   const { return loc_topol_cmpl_;  }
  const PS::ReallocatableArray<PS::U32>& loc_topol_imcmpl() const { return loc_topol_imcmpl_; }

  explicit ForceBondedMPI(const PS::U32 est_loc_amp) {
    ampid.resizeNoInitialize(est_loc_amp);
    ampid_buf.resizeNoInitialize(est_loc_amp);
    
    loc_topol_cmpl_.resizeNoInitialize(est_loc_amp * Parameter::all_unit);
    loc_topol_imcmpl_.resizeNoInitialize(est_loc_amp * Parameter::all_unit);
    is_real_surf.resizeNoInitialize(est_loc_amp * Parameter::all_unit);

    cf_bend_cmpl.resizeNoInitialize(est_loc_amp);
    cf_bend_imcmpl.resizeNoInitialize(est_loc_amp);
  }
  ~ForceBondedMPI() {}

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
    
    if (mask[0]) {
      F[0] -= Fbond;
      d_vir.x += dr.x * Fbond.x;
      d_vir.y += dr.y * Fbond.y;
      d_vir.z += dr.z * Fbond.z;
    }
    if (mask[1]) {
      F[1] += Fbond;
      d_vir.x += dr.x * Fbond.x;
      d_vir.y += dr.y * Fbond.y;
      d_vir.z += dr.z * Fbond.z;
    }
  }

  static void StoreBendForceNoARLaw(const PS::F64vec*	__restrict dr,
				    const PS::F64*	__restrict inv_dr,
				    const PS::F64*	__restrict dist,
				    PS::F64vec&		__restrict d_vir,
				    PS::F64&		__restrict d_lap,
				    PS::F64vec*		__restrict F,
				    const bool*		__restrict mask,
				    const PS::F64       cf_b)
  {
    const PS::F64	inv_dr_prod	= inv_dr[0] * inv_dr[1];
    const PS::F64	inv_dist[2]	= { inv_dr[0] * inv_dr[0],
					    inv_dr[1] * inv_dr[1] };
    const PS::F64	in_prod		= dr[0] * dr[1];
    const PS::F64	cf_bd		= cf_b * inv_dr_prod;
    const PS::F64       cf_crs[2]	= { in_prod * inv_dist[0],
					    in_prod * inv_dist[1] };
    
    const PS::F64vec Ftb0(cf_bd * (dr[1].x - cf_crs[0] * dr[0].x),
			  cf_bd * (dr[1].y - cf_crs[0] * dr[0].y),
			  cf_bd * (dr[1].z - cf_crs[0] * dr[0].z));
    const PS::F64vec Ftb1(cf_bd * (dr[0].x - cf_crs[1] * dr[1].x),
			  cf_bd * (dr[0].y - cf_crs[1] * dr[1].y),
			  cf_bd * (dr[0].z - cf_crs[1] * dr[1].z));

    if (mask[0]) {
      F[0] -= Ftb0;
      d_vir.x += dr[0].x * Ftb0.x;
      d_vir.y += dr[0].y * Ftb0.y;
      d_vir.z += dr[0].z * Ftb0.z;
    }
    if (mask[1]) {
      F[1] += Ftb0 - Ftb1;
      d_vir.x += dr[0].x * Ftb0.x + dr[1].x * Ftb1.x;
      d_vir.y += dr[0].y * Ftb0.y + dr[1].y * Ftb1.y;
      d_vir.z += dr[0].z * Ftb0.z + dr[1].z * Ftb1.z;
    }
    if (mask[2]) {
      F[2] += Ftb1;
      d_vir.x += dr[1].x * Ftb1.x;
      d_vir.y += dr[1].y * Ftb1.y;
      d_vir.z += dr[1].z * Ftb1.z;
    }
  }

  //NOTE: minimum image convention is not needed.
  //      this function is designed for small # of bonds.
  //ASSUME: bond_n >= 3.
  template<PS::U32 bond_n>
  void CalcBondBendLocalCell(Tpsys&		__restrict sys,
			     const PS::U32	beg_bond_id,
			     const PS::F64      cf_b,
			     PS::F64vec&	__restrict d_vir,
			     PS::F64&		__restrict d_lap)
  {
    PS::F64vec Fbb[bond_n], pos_buf[bond_n], dr[bond_n - 1];
    PS::F64  dist2[bond_n - 1], inv_dr[bond_n - 1];

    const PS::U32* l_dst = &(loc_topol_cmpl_[beg_bond_id]);
    
    pos_buf[0] = sys[ l_dst[0] ].pos;
    pos_buf[1] = sys[ l_dst[1] ].pos;

    dr[0] = pos_buf[1] - pos_buf[0];
    ForceBonded<Tpsys>::MinImage(dr[0]);
    dist2[0] = dr[0] * dr[0];
    inv_dr[0] = 1.0 / std::sqrt(dist2[0]);
    
    ForceBonded<Tpsys>::StoreBondForceWithARLaw(dr[0], inv_dr[0], d_vir, d_lap, &Fbb[0]);
#pragma unroll
    for (PS::U32 unit = 2; unit < bond_n; unit++) {
      pos_buf[unit] = sys[ l_dst[unit] ].pos;
      dr[unit - 1] = pos_buf[unit] - pos_buf[unit - 1];
      ForceBonded<Tpsys>::MinImage(dr[unit - 1]);
      dist2[unit - 1] = dr[unit - 1] * dr[unit - 1];
      inv_dr[unit - 1] = 1.0 / std::sqrt(dist2[unit - 1]);
      
      ForceBonded<Tpsys>::StoreBondForceWithARLaw(dr[unit - 1], inv_dr[unit - 1], d_vir, d_lap, &Fbb[unit - 1]);
      ForceBonded<Tpsys>::StoreBendForceWithARLaw(&dr[unit - 2], &inv_dr[unit - 2], dist2, d_vir, d_lap, &Fbb[unit - 2], cf_b);
    }

    //Store the sum of force.
#pragma unroll
    for (PS::U32 unit = 0; unit < bond_n; unit++)
      sys[ l_dst[unit] ].acc += Fbb[unit];
  }

  template<PS::U32 bond_n>
  void CalcBondBendSurface(Tpsys& __restrict sys,
			   const PS::ReallocatableArray<Pepj>& epj_org,
			   const PS::U32 beg_bond_id,
			   const PS::F64 cf_b,
			   PS::F64vec& __restrict d_vir,
			   PS::F64&    __restrict d_lap)
  {
    PS::F64vec Fbb[bond_n], pos_buf[bond_n], dr[bond_n - 1];
    PS::F64  dist2[bond_n - 1], inv_dr[bond_n - 1];

    const bool* mask = &(is_real_surf[beg_bond_id]);
    const PS::U32* l_dst = &(loc_topol_imcmpl_[beg_bond_id]);
    
    if (l_dst[0] != 0xffffffff) pos_buf[0] = epj_org[ l_dst[0] ].pos;
    if (l_dst[1] != 0xffffffff) pos_buf[1] = epj_org[ l_dst[1] ].pos;

    dr[0] = pos_buf[1] - pos_buf[0];
    ForceBonded<Tpsys>::MinImage(dr[0]);
    dist2[0] = dr[0] * dr[0];
    inv_dr[0] = 1.0 / std::sqrt(dist2[0]);
    
    StoreBondForceNoARLaw(dr[0], inv_dr[0], d_vir, d_lap, &Fbb[0], &mask[0]);
#pragma unroll
    for (PS::U32 unit = 2; unit < bond_n; unit++) {
      if (l_dst[unit] != 0xffffffff) pos_buf[unit] = epj_org[ l_dst[unit] ].pos;
      dr[unit - 1] = pos_buf[unit] - pos_buf[unit - 1];
      ForceBonded<Tpsys>::MinImage(dr[unit - 1]);
      dist2[unit - 1] = dr[unit - 1] * dr[unit - 1];
      inv_dr[unit - 1] = 1.0 / std::sqrt(dist2[unit - 1]);

      StoreBondForceNoARLaw(dr[unit - 1], inv_dr[unit - 1], d_vir, d_lap, &Fbb[unit - 1], &mask[unit - 1]);
      StoreBendForceNoARLaw(&dr[unit - 2], &inv_dr[unit - 2], dist2, d_vir, d_lap, &Fbb[unit - 2], &mask[unit - 2], cf_b);
    }

    //Store the sum of force.
#pragma unroll
    for (PS::U32 unit = 0; unit < bond_n; unit++)
      if (mask[unit]) sys[ l_dst[unit] ].acc += Fbb[unit];
  }

  void AppendCfBend(PS::ReallocatableArray<PS::F64>& cf_bend,
		    const std::vector<PS::U32>& core_amp_id,
		    const PS::U32 cur_id) {
    PS::F64 cf_temp = Parameter::cf_b;
    for (PS::U32 j = 0; j < core_amp_id.size(); j++) {
      if (cur_id == core_amp_id[j]) cf_temp = Parameter::cf_b_rigid;
    }
    cf_bend.pushBackNoCheck(cf_temp);
  }

  void MakeLocalBondedList(const Tpsys& sys,
			   const PS::ReallocatableArray<Pepj>& epj_org,
			   const std::vector<PS::U32>& core_amp_id) {
    const PS::U32 real_n = sys.getNumberOfParticleLocal();
    const PS::U32 all_n  = epj_org.size();

    // resize buffer if needed.
    ampid.resizeNoInitialize(all_n);
    ampid_buf.resizeNoInitialize(all_n);
    loc_topol_cmpl_.resizeNoInitialize(all_n);
    loc_topol_imcmpl_.resizeNoInitialize(all_n * Parameter::all_unit);
    is_real_surf.resizeNoInitialize(10 * all_n * Parameter::all_unit);
    cf_bend_cmpl.resizeNoInitialize(all_n);
    cf_bend_imcmpl.resizeNoInitialize(all_n);
    
    // copy to buffer.
    ampid.clearSize();
    ampid2idx buf;
    for (PS::U32 i = 0; i < all_n; i++) {
      if (epj_org[i].prop != Parameter::Solvent) {
	buf.amp_id	= epj_org[i].amp_id;
	buf.unit	= epj_org[i].unit;
	buf.idx		= i;
	buf.is_real	= (i < real_n);
	buf.setKey();
	ampid.pushBackNoCheck(buf);
      }
    }

    // remove duplications.
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
    using __gnu_parallel::sort;
#else
    using std::sort;
#endif
    PS::U32 amp_ptcl_loc = ampid.size();
    sort(ampid.getPointer(),
	 ampid.getPointer() + amp_ptcl_loc,
	 [](const ampid2idx& i, const ampid2idx& j) {
	   const PS::U32 ikey = i.key, jkey = j.key;
	   if (ikey != jkey)
	     return (ikey < jkey);
	   else
	     return (i.is_real > j.is_real);
	 });
    const auto new_end = std::unique(ampid.getPointer(),
				     ampid.getPointer() + amp_ptcl_loc,
				     [](const ampid2idx& i, const ampid2idx& j) {
				       return i.key == j.key;
				     });
    amp_ptcl_loc = std::distance(ampid.getPointer(), new_end);
    ampid[amp_ptcl_loc].amp_id = 0xffffffff; // add dummy id

    // then construct complete and imcomplete topology.
    cmplt_amp_ = imcmplt_amp_ = 0;
    cf_bend_cmpl.clearSize(); cf_bend_imcmpl.clearSize();
    PS::U32 id_cmpl = 0, id_imcmpl = 0, cnt_real = ampid[0].is_real, cnt = 1, id_bef = ampid[0].amp_id, id_cur;
    PS::U32 imcmpl_buf[Parameter::all_unit] = {0xffffffff};
    bool    isreal_buf[Parameter::all_unit] = {false};
    for (PS::U32 i = 1; i < amp_ptcl_loc + 1; i++) {
      id_cur = ampid[i].amp_id;

      if (id_cur == id_bef) {
	cnt_real += ampid[i].is_real;
	cnt++;
      } else {
	if(cnt_real == Parameter::all_unit) {
	  const PS::U32 beg = i - Parameter::all_unit;
	  for (PS::U32 j = 0; j < Parameter::all_unit; j++) {
	    loc_topol_cmpl_[id_cmpl++] = ampid[beg + j].idx;
	  }
	  AppendCfBend(cf_bend_cmpl, core_amp_id, id_bef);
	  cmplt_amp_++;
	} else {
	  //initialize buffer
	  for (PS::U32 j = 0; j < Parameter::all_unit; j++) {
	    imcmpl_buf[j] = 0xffffffff;
	    isreal_buf[j] = false;
	  }
	  
	  const PS::U32 beg = i - cnt;
	  for (PS::U32 j = 0; j < cnt; j++) {
	    const PS::U32 dst	= ampid[beg + j].unit;
	    imcmpl_buf[dst]	= ampid[beg + j].idx;
	    isreal_buf[dst]	= ampid[beg + j].is_real;
	  }

	  for (PS::U32 j = 0; j < Parameter::all_unit; j++) {
	    loc_topol_imcmpl_[id_imcmpl] = imcmpl_buf[j];
	    is_real_surf[id_imcmpl]	= isreal_buf[j];
	    id_imcmpl++;
	  }

	  AppendCfBend(cf_bend_imcmpl, core_amp_id, id_bef);
	  imcmplt_amp_++;
	}
	cnt_real = (ampid[i].is_real);
	cnt = 1;
      }
      id_bef = id_cur;
    }
  }

  void CheckSurfaceTopol() const {
    PS::U32 cnt = 0;
    for (PS::U32 aid = 0; aid < imcmplt_amp_; aid++) {
      bool req_flag[Parameter::all_unit] = { false };
      for (PS::U32 unit = 0; unit < Parameter::all_unit; unit++) {
	if (is_real_surf[cnt]) {
	  for(PS::S32 j = -2; j < 3; j++) {
	    const PS::S32 unit_j = unit + j;
	    if(unit_j >= 0 && unit_j < static_cast<PS::S32>(Parameter::all_unit) )
	      req_flag[unit_j] = true;
	  }
	}
	cnt++;
      }
      
      cnt -= Parameter::all_unit;

      for (PS::U32 unit = 0; unit < Parameter::all_unit; unit++) {
	if (req_flag[unit] && (loc_topol_imcmpl_[cnt] == 0xffffffff)) {
	  std::cerr << "Missing pair particles\n";
	  PS::Abort();
	}
	cnt++;
      }
    }
  }

  void CalcListedForce(Tpsys& sys,
		       const PS::ReallocatableArray<Pepj>& epj_org,
		       PS::F64vec& bonded_vir,
		       const std::vector<PS::U32>& core_amp_id) {
    MakeLocalBondedList(sys, epj_org, core_amp_id);
    CheckSurfaceTopol();
    
    PS::F64vec d_vir(0.0); PS::F64 d_lap = 0.0;
    for (PS::U32 i = 0; i < cmplt_amp_; i++) {
      CalcBondBendLocalCell<Parameter::all_unit>(sys, i * Parameter::all_unit, cf_bend_cmpl[i], d_vir, d_lap);
    }

    //for inter cell
    for (PS::U32 i = 0; i < imcmplt_amp_; i++) {
      CalcBondBendSurface<Parameter::all_unit>(sys, epj_org, i * Parameter::all_unit, cf_bend_imcmpl[i], d_vir, d_lap);
    }

    bonded_vir = d_vir;
  }

  void CheckCompleteTopol(const Tpsys& sys) const {
    for (PS::U32 i = 0; i < cmplt_amp_; i++) {
      const PS::U32 amp_id = sys[loc_topol_cmpl_[Parameter::all_unit * i]].amp_id;
      for (PS::U32 j = 1; j < Parameter::all_unit; j++) {
	assert(loc_topol_cmpl_[Parameter::all_unit * i + j] < static_cast<PS::U32>(sys.getNumberOfParticleLocal()));
	if (amp_id != sys[ loc_topol_cmpl_[Parameter::all_unit * i + j] ].amp_id) {
	  std::cerr << "Fail topology build\n";
	  std::cerr << amp_id << " " << sys[ loc_topol_cmpl_[Parameter::all_unit * i + j] ].amp_id << std::endl;
	  std::cerr << loc_topol_cmpl_[Parameter::all_unit * i + j] << " " << sys.getNumberOfParticleLocal() << std::endl;
	  std::cerr << i << " " << j << std::endl;
	} else {
	  assert(loc_topol_cmpl_[Parameter::all_unit * i + j] < static_cast<PS::U32>(sys.getNumberOfParticleLocal()));
	}
      }	
    }
  }
};
