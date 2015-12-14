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
  
};

struct FPDPD {
  PS::S64 id;
  PS::S32 prop, unit, time;
  PS::F64vec pos;
  PS::F64vec vel, vel_bef;
  PS::F64vec acc, acc_bef;
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
    
  }
  void writeAscii(FILE *fp) {
    
  }
  
  
};

struct EPIDPD {
  PS::S64 id;
  PS::S32 prop, time;
  PS::F64vec pos, vel;

  PS::F64vec getPos() const {
    return this->pos;
  }

  void copyFromFP(const FPDPD& fp) {
    this->pos = fp.pos;
    this->vel = fp.vel;
    this->id = fp.id;
    this->prop = fp.prop;
    this->time = fp.time;
  }
};

struct EPJDPD {

};

#define SARU(ix,iy,iz) Saru saru( (ix) , (iy) , (iz) )
#define CALL_SARU(x, y) \
  std::sqrt( -2.0f * std::log(saru.f( (x), (y) ) ) ) * std::cos(2.0f * M_PI * saru.f( (x), (y)) )

struct CalcForceEpEp {
  //interactions
  static PS::F64 cf_c[Parameter::prop_num][Parameter::prop_num];
  static PS::F64 cf_g[Parameter::prop_num][Parameter::prop_num];
  static PS::F64 cf_r[Parameter::prop_num][Parameter::prop_num];

  //cutoff length
  static PS::F64 rc, rc2;

  //prng seed
  static PS::S64 m_seed;
  
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
      const PS::S64   idi = epi[i].id;
      const PS::S32 propi = epi[i].prop;
      
      PS::F64vec fsum(0.0, 0.0, 0.0);
      PS::F64vec psum(0.0, 0.0, 0.0);
      
      for(PS::S32 j = 0; j < nj; j++) {
	if(idi == epj[j].id) continue;
	
	PS::F64vec drij = ri - epj[j].pos;
	PS::F64 dr2 = drij * drij;
	if(dr2 < rc2) {
	  const PS::F64vec dvij = vi - epj[j].vel;
	  
	  const PS::S64 idj = epj[j].id;
	  const PS::S32 time = epi[i].time;
	  
	  PS::S64 m_i = idi, m_j = idj;
	  if(idi > idj){ // m_i <= m_j 
	    m_i = idj;
	    m_j = idi;
	  }
	  
	  SARU(m_i, m_j, m_seed + time);
	  const PS::F64 nrml_rnd = CALL_SARU(-1, 1);
	  
	  const PS::F64 dr = std::sqrt(dr2);
	  const PS::F64 inv_dr = 1.0 / dr;
	  const PS::S32 propj  = epj[j].prop;
	  const PS::F64 one_m_dr = 1.0 - drij * inv_dr;

	  const PS::F64 wrij = one_m_dr; // pow = 1
	  //const PS::F64 wrij = std::sqrt(one_m_dr); //pow = 1 / 2
	  const PS::F64 sq_wrij = std::sqrt(one_m_dr);
	  
	  const PS::F64 drij_dvij = drij * dvij;
	  
	  const PS::F64 all_cf = ( cf_c[propi][propj] * one_m_dr - cf_g[propi][propj] * wrij * drij_dvij * inv_dr 
				   + cf_r[propi][propj] * sq_wrij * nrml_rnd ) * inv_dr;

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

