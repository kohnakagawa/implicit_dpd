#pragma once

#include "saruprng.hpp"

namespace RESULT {
  struct ForceDPD {
    PS::F64vec acc;
    PS::F64vec press;
  
    void clear(){
      acc   = 0.0;
      press = 0.0;
    }
  };

  struct Density {
    PS::F64 dens[2];
    void clear() {
      dens[0] = 0.0;
      dens[1] = 0.0;
    }
  };
};

struct FPDPD {
  PS::U32 id, prop, amp_id, unit;
  PS::F64vec pos;
  PS::F64vec vel, vel_buf;
  PS::F64vec acc;
  PS::F64vec press;
  PS::F64 density[2];

  //essential member functions
  PS::F64 getRSearch() const {
    return Parameter::search_rad;
  }
  PS::F64vec getPos() const {
    return this->pos;
  }
  void setPos(const PS::F64vec& p) {
    pos = p;
  }
  void copyFromForce(const RESULT::ForceDPD& force) {
    acc = force.acc;
    press = force.press;
  }

  void copyFromForce(const RESULT::Density& dens) {
    density[0] = dens.dens[0];
    density[1] = dens.dens[1];
  }
  
  //for I/O
  void readAscii(FILE *fp) {
    fscanf(fp, "%u %u %u %u %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	   &id, &prop, &amp_id, &unit,
	   &(pos.x), &(pos.y), &(pos.z),
	   &(vel.x), &(vel.y), &(vel.z), &(vel_buf.x), &(vel_buf.y), &(vel_buf.z),
	   &(acc.x), &(acc.y), &(acc.z),
	   &(density[0]), &(density[1]));
  }
  void writeAscii(FILE *fp) const {
    fprintf(fp, "%u %u %u %u %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n",
	    id, prop, amp_id, unit,
	    pos.x, pos.y, pos.z,
	    vel.x, vel.y, vel.z, vel_buf.x, vel_buf.y, vel_buf.z,
	    acc.x, acc.y, acc.z,
	    density[0], density[1]);
  }
};

namespace EPI {
  struct DPD {
    PS::U32 id, prop;
    PS::F64vec pos, vel;
    PS::F64 dens[2];

    PS::F64vec getPos() const {
      return this->pos;
    }

    PS::F64 getRSearch() const {
      return Parameter::search_rad;
    }

    void copyFromFP(const FPDPD& fp) {
      this->pos		= fp.pos;
      this->vel		= fp.vel;
      this->id		= fp.id;
      this->prop	= fp.prop;
      this->dens[0]	= fp.density[0];
      this->dens[1]	= fp.density[1];
    }

  };
  
  struct Density {
    PS::U32 prop;
    PS::F64vec pos;

    PS::F64vec getPos() const {
      return this->pos;
    }

    PS::F64 getRSearch() const {
      return Parameter::search_rad;
    }
    
    void copyFromFP(const FPDPD& fp) {
      this->prop	= fp.prop;
      this->pos		= fp.pos;
    }
  };
};

namespace EPJ {
  struct DPD {
    PS::U32 id, prop;
    PS::F64vec pos, vel;
    PS::F64 dens[2];
    
    void copyFromFP(const FPDPD& fp) {
      this->id			= fp.id;
      this->prop		= fp.prop;
      this->pos			= fp.pos;
      this->vel			= fp.vel;
      this->dens[0]		= fp.density[0];
      this->dens[1]		= fp.density[1];
    }
    
    PS::F64vec getPos() const {
      return this->pos;
    }
    
    void setPos(const PS::F64vec& pos_) {
      this->pos = pos_;
    }
  };

  struct Density {
    PS::U32 prop;
    PS::F64vec pos;
    
    void copyFromFP(const FPDPD& fp) {
      this->prop		= fp.prop;
      this->pos			= fp.pos;
    }
    
    PS::F64vec getPos() const {
      return this->pos;
    }
    
    void setPos(const PS::F64vec& pos_) {
      this->pos = pos_;
    }
  };
};

#define SARU(ix,iy,iz) Saru saru( (ix) , (iy) , (iz) )
#define CALL_SARU_NRML(x, y) \
  std::sqrt( -2.0f * std::log(saru.f( (x), (y) ) ) ) * std::cos(2.0f * M_PI * saru.f( (x), (y)) )
#define CALL_SARU(x, y) \
  saru.f( (x), (y) )

struct CalcDensity {
  void operator () (const EPI::Density* __restrict epi,
		    const PS::S32 ni,
		    const EPJ::Density* __restrict epj,
		    const PS::S32 nj,
		    RESULT::Density* __restrict result)
  {
    for(PS::S32 i = 0; i < ni; i++) {
      const PS::F64vec ri = epi[i].pos;
      PS::F64 d_sum[2] = { 0.0 };
      for(PS::S32 j = 0; j < nj; j++) {
	const PS::U32 propj = epj[i].prop;
	
	const PS::F64vec drij = ri - epj[j].pos;
	const PS::F64 dr2 = drij * drij;
	if(dr2 < Parameter::rc2 && propj < 2) {
	  const PS::F64 dr = std::sqrt(dr2);
	  d_sum[propj] += (Parameter::rc - dr) * (Parameter::rc - dr); //NOTE: density kernel is harmonic
	}
      }
      result[i].dens[0] += d_sum[0];
      result[i].dens[1] += d_sum[1];
    }
  }
};

struct CalcForceEpEpDPD {
  //prng seed
  static PS::U32 m_seed;
  
  //NOTE: Is the minimum image convention employed?
  //ASSUME: particles are hydrophilic or hydrophobic.
  void operator () (const EPI::DPD* __restrict epi,
		    const PS::S32 ni,
		    const EPJ::DPD* __restrict epj,
		    const PS::S32 nj,
		    RESULT::ForceDPD* __restrict result)
  {
    for(PS::S32 i = 0; i < ni; i++) {
      const PS::F64vec ri = epi[i].pos;
      const PS::F64vec vi = epi[i].vel;
      const PS::U32   idi = epi[i].id;
      const PS::U32 propi = epi[i].prop;
      const PS::F64 densi[] = {epi[i].dens[0], epi[i].dens[1]};
      
      PS::F64vec fsum(0.0, 0.0, 0.0);
      PS::F64vec psum(0.0, 0.0, 0.0);
      
      for(PS::S32 j = 0; j < nj; j++) {
	if(idi == epj[j].id) continue;
	
	const PS::F64vec drij = ri - epj[j].pos;
	const PS::F64 dr2 = drij * drij;
	if(dr2 < Parameter::rc2) {
	  const PS::F64vec dvij = vi - epj[j].vel;
	  const PS::F64  densj[]  = {epj[j].dens[0], epj[j].dens[1]};
	  
	  const PS::U32 idj = epj[j].id;
	  
	  PS::U32 m_i = idi, m_j = idj;
	  if(idi > idj){ // m_i <= m_j 
	    m_i = idj;
	    m_j = idi;
	  }
	  
	  SARU(m_i, m_j, m_seed + Parameter::time);
	  const PS::F64 rnd = CALL_SARU(-1, 1); //uniform 
	  //const PS::F64 rnd = CALL_SARU_NRML(-1, 1); //normal

	  //kernel for thermostat
	  const PS::F64 dr = std::sqrt(dr2);
	  const PS::F64 inv_dr = 1.0 / dr;
	  const PS::U32 propj  = epj[j].prop;
	  const PS::F64 one_m_dr = 1.0 - dr * Parameter::irc;

	  //kernel for conservative force
	  const PS::F64 cf_co  = Parameter::cf_c[propi][propj] * 6.0 * (dr - Parameter::arc) * (dr - Parameter::rc) * (dr >= Parameter::arc);
	  const PS::F64 cf_mbd = ( (densi[0] + densj[0]) * Parameter::cf_m[propi][propj][0] +
				   (densi[1] + densj[1]) * Parameter::cf_m[propi][propj][1] ) * one_m_dr;

	  const PS::F64 wrij = one_m_dr; // pow = 1
	  //const PS::F64 wrij = std::sqrt(one_m_dr); //pow = 1 / 2
	  const PS::F64 sq_wrij = std::sqrt(wrij);
	  
	  const PS::F64 drij_dvij = drij * dvij;
	  
	  const PS::F64 all_cf = (cf_co + cf_mbd +  //conservative
				  Parameter::cf_r[propi][propj] * sq_wrij * rnd - //random
				  Parameter::cf_g[propi][propj] * wrij * drij_dvij * inv_dr) * inv_dr; //dissipation

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

PS::U32 CalcForceEpEpDPD::m_seed = 0;
