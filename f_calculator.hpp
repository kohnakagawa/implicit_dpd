#pragma once

#include "../../src/particle_simulator.hpp"
#include "saruprng.hpp"
#include "parameter.hpp"

struct ForceDPD {
  PS::F64vec acc;
  PS::F64vec press;
  
  void clear(){
    acc   = 0.0;
    press = 0.0;
  }
};

struct ForceBonded {

  struct ampid2idx {
    PS::U32 amp_id, unit, idx;
    inline PS::U32 getKey() { return Parameter::all_unit * amp_id + unit; }
  };

  PS::U32 cmplt_amp = 0, imcmplt_amp = 0;
  PS::ReallocatableArray<ampid2idx> ampid, ampid_buf;
  PS::ReallocatableArray<PS::U32> loc_topol_cmpl, loc_topol_imcmpl;
  PS::RadixSort<PS::U32> rsorter;

  ForceBonded(const int est_loc_amp) {
    ampid.resizeNoInitialize(est_loc_amp);
    ampid_buf.resizeNoInitialize(est_loc_amp);
    
    loc_topol_cmpl.resizeNoInitialize(est_loc_amp * Parameter::all_unit);
    loc_topol_imcmpl.resizeNoInitialize(est_loc_amp * Parameter::all_unit);
  }
  ~ForceBonded() {
    
  }

  void StoreBondForceNoARLaw( ) {
    
  }

  void StoreBondForceWithARLaw(const PS::F64vec& __restrict dr,
			       const PS::F64&    __restrict inv_dr,
			       PS::F64vec& __restrict d_virial,
			       PS::F64& __restrict lap_conf,
			       PS::F64vec* __restrict F)
  {
    const double cf_bond = Parameter::cf_s * (inv_dr - Parameter::ibond);
    
    const PS::F64vec Fbond(cf_bond * dr.x, cf_bond * dr.y, cf_bond * dr.z);
    
    d_virial.x += dr.x * Fbond.x;
    d_virial.y += dr.y * Fbond.y;
    d_virial.z += dr.z * Fbond.z;

    lap_conf += Parameter::cf_s * (6.0 * Parameter::ibond - 4.0 * inv_dr);
    
    F[0] -= Fbond;
    F[1] += Fbond;
  }

  void StoreBendForceWithARLaw(const PS::F64vec* __restrict dr,
			       const PS::F64* __restrict inv_dr,
			       const PS::F64* __restrict dist,
			       PS::F64vec& __restrict d_virial,
			       PS::F64& __restrict lap_conf,
			       PS::F64vec* __restrict F)
  {
    const PS::F64	inv_dr_prod	= inv_dr[0] * inv_dr[1];
    const PS::F64	inv_dist[2]	= {inv_dr[0] * inv_dr[0],
					   inv_dr[1] * inv_dr[1] };
    const PS::F64	in_prod		= dr[0] * dr[1];
    const PS::F64	cf_b		= itrs.cf_bend * inv_dr_prod;
    const PS::F64       cf_crs[2]	= {in_prod * inv_dist[0],
					   in_prod * inv_dist[1] };
    
    const PS::F64vec Ftb0(cf_b * (dr[1].x - cf_crs[0]*dr[0].x),
			  cf_b * (dr[1].y - cf_crs[0]*dr[0].y),
			  cf_b * (dr[1].z - cf_crs[0]*dr[0].z));
    const PS::F64vec Ftb1(cf_b * (dr[0].x - cf_crs[1]*dr[1].x),
			  cf_b * (dr[0].y - cf_crs[1]*dr[1].y),
			  cf_b * (dr[0].z - cf_crs[1]*dr[1].z));
    
    d_virial.x += dr[0].x * Ftb0.x + dr[1].x * Ftb1.x;
    d_virial.y += dr[0].y * Ftb0.y + dr[1].y * Ftb1.y;
    d_virial.z += dr[0].z * Ftb0.z + dr[1].z * Ftb1.z;

    lap_conf += 2.0 * cf_b * inv_dist[0] * inv_dist[1] * ( in_prod * ( in_prod + 2.0 * (dist[0] + dist[1]) ) + dist[0] * dist[1]);
    
    F[0] -= Ftb0;
    F[1] += Ftb0 - Ftb1;
    F[2] += Ftb1;
  }

  template<class Tpsys, int bond_n>
  void CalcBondBendNoConfigTemp(Tpsys& sys) {
    PS::F64vec Fbb[bond_n], temp_pos[bond_n], dr[bond_n - 1];
    PS::F64 dist2[bond_n - 1], inv_dr[bond_n - 1];

    const int l_dst[2] = {elem_idx[beg_idx],elem_idx[beg_idx+1]};
    temp_pos[0] = pr[l_dst[0]]; temp_pos[1] = pr[l_dst[1]];
    dr[0] = temp_pos[1] - temp_pos[0];
    MinImage(dr[0],param);
    dist2[0] = dr[0].dist2();
    inv_dr[0] = 1.0 / sqrt(dist2[0]);
    
    StoreBondForce(dr[0], inv_dr[0], d_virial, lap_conf, &Fbb[0], param);

    for(int unit=2; unit < bond_n; unit++){
      const int load_dest = elem_idx[beg_idx + unit];
      temp_pos[unit] = pr[load_dest];
      dr[unit-1] = temp_pos[unit] - temp_pos[unit-1];
      dist2[unit-1] = dr[unit-1].dist2();
      inv_dr[unit-1] = 1.0 / std::sqrt(dist2[unit-1]);
      
      StoreBondForce(dr[unit-1], inv_dr[unit-1], d_virial, lap_conf, &Fbb[unit-1], param);
      StoreBendForce(&dr[unit-2], &inv_dr[unit-2], dist2, d_virial, lap_conf, &Fbb[unit-2], param);
    }

    for(int unit = 0; unit < bond_n; unit++) {
      const int str_dest = elem_idx[beg_idx + unit];
      force[str_dest] += Fbb[unit];
    }
  }
  
  
  template<class Tpsys>
  void MakeLocalBondedList(const Tpsys& sys) {
    PS::U32 n = sys.getNumberOfParticleLocal();
    
    for(PS::U32 i = 0; i < n; i++) {
      ampid[i].amp_id	= sys[i].amp_id;
      ampid[i].unit	= sys[i].unit;
      ampid[i].idx	= i;
    }
    rsorter.lsdSort(ampid, ampid_buf, 0, n);
    
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
	}
	cnt = 1;
      }
      id_bef = id_cur;
    }
  }
  
  template<class Tpsys>
  void CalcListedForce(Tpsys& sys) {
    for(PS::U32 i = 0; i < cmplt_amp; i++) {
      
      
    }

    //
    for(PS::U32 i = 0; i < imcmplt_amp; i++) {
      
    }
    
  }
  
};

struct FPDPD {
  PS::U32 id, prop, amp_id, unit;
  PS::F64vec pos;
  PS::F64vec vel, vel_buf;
  PS::F64vec acc;
  PS::F64vec press;
  PS::F64 search_radius;
  
  //essential member functions
  PS::F64vec getPos() const {
    return this->pos;
  }
  void copyFromForce(const ForceDPD& force) {
    acc = force.acc;
    press = force.press;
  }
  PS::F64 getRSearch() const {
    return this->search_radius;
  }
  
  //for I/O
  void readAscii(FILE *fp) {
    fprintf(fp, "%u %u %u %u %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n",
	    &id, &prop, &amp_id, &unit,
	    &pos.x, &pos.y, &pos.z,
	    &vel.x, &vel.y, &vel.z, &vel_buf.x, &vel_buf.y, &vel_buf.z,
	    &acc.x, &acc.y, &acc.z,);
  }
  void writeAscii(FILE *fp) const {
    fprintf(fp, "%u %u %u %u %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n",
	    id, prop, amp_id, unit,
	    pos.x, pos.y, pos.z,
	    vel.x, vel.y, vel.z, vel_buf.x, vel_buf.y, vel_buf.z,
	    acc.x, acc.y, acc.z);
  }
  
  
};

struct EPIDPD {
  PS::U32 id, prop;
  PS::F64vec pos, vel;

  PS::F64vec getPos() const {
    return this->pos;
  }

  void copyFromFP(const FPDPD& fp) {
    this->pos = fp.pos;
    this->vel = fp.vel;
    this->id  = fp.id;
    this->prop = fp.prop;
  }
};

struct EPJDPD {

};

#define SARU(ix,iy,iz) Saru saru( (ix) , (iy) , (iz) )
#define CALL_SARU_NRML(x, y) \
  std::sqrt( -2.0f * std::log(saru.f( (x), (y) ) ) ) * std::cos(2.0f * M_PI * saru.f( (x), (y)) )
#define CALL_SARU(x, y) \
  saru.f( (x), (y) )

struct CalcForceEpEp {
  //prng seed
  static PS::U32 m_seed;
  
  //NOTE: Is the minimum image convention employed?
  void operator () (const EPIDPD *epi,
		    const PS::S32 ni,
		    const EPJDPD *epj,
		    const PS::S32 nj,
		    ForceDPD * result)
  {
    for(PS::S32 i = 0; i < ni; i++) {
      const PS::F64vec ri = epi[i].pos;
      const PS::F64vec vi = epi[i].vel;
      const PS::U32   idi = epi[i].id;
      const PS::U32 propi = epi[i].prop;
      
      PS::F64vec fsum(0.0, 0.0, 0.0);
      PS::F64vec psum(0.0, 0.0, 0.0);
      
      for(PS::S32 j = 0; j < nj; j++) {
	if(idi == epj[j].id) continue;
	
	const PS::F64vec drij = ri - epj[j].pos;
	const PS::F64 dr2 = drij * drij;
	if(dr2 < Parameter::rc2) {
	  const PS::F64vec dvij = vi - epj[j].vel;
	  
	  const PS::U32 idj = epj[j].id;
	  
	  PS::U32 m_i = idi, m_j = idj;
	  if(idi > idj){ // m_i <= m_j 
	    m_i = idj;
	    m_j = idi;
	  }
	  
	  SARU(m_i, m_j, m_seed + Parameter::time);
	  const PS::F64 rnd = CALL_SARU(-1, 1); //uniform 
	  //const PS::F64 rnd = CALL_SARU_NRML(-1, 1); //normal
	  
	  const PS::F64 dr = std::sqrt(dr2);
	  const PS::F64 inv_dr = 1.0 / dr;
	  const PS::U32 propj  = epj[j].prop;
	  const PS::F64 one_m_dr = 1.0 - drij * inv_dr;

	  const PS::F64 wrij = one_m_dr; // pow = 1
	  //const PS::F64 wrij = std::sqrt(one_m_dr); //pow = 1 / 2
	  const PS::F64 sq_wrij = std::sqrt(one_m_dr);
	  
	  const PS::F64 drij_dvij = drij * dvij;
	  
	  const PS::F64 all_cf = ( Parameter::cf_c[propi][propj] * one_m_dr - Parameter::cf_g[propi][propj] * wrij * drij_dvij * inv_dr 
				   + Parameter::cf_r[propi][propj] * sq_wrij * rnd ) * inv_dr;

	  const PS::F64vec dF(all_cf * drij.x, all_cf * drij.y, all_cf * drij.z);
	  const PS::F64vec dP(dF.x * drij.x, dF.y * drij.y, dF.z * drij.z);
	  
	  fsum += dF;
	  psum += dP;
	}
      }
      result[i].acc += fsum;
      result[i].press += psum;
    }
  }
};

#undef SARU
#undef CALL_SARU_NRML
#undef CALL_SARU
