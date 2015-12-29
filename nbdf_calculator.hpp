#pragma once

#include "saruprng.hpp"
// #include "static_for.hpp"
#include "ptcl_class.hpp"

struct CalcDensity {
  void operator () (const EPI::Density* epi,
		    const PS::S32 ni,
		    const EPJ::Density* epj,
		    const PS::S32 nj,
		    RESULT::Density* result)
  {
    for(PS::S32 i = 0; i < ni; i++) {
      const PS::F64vec ri = epi[i].pos;
      PS::F64 d_sum[Parameter::prop_num] = { 0.0 };
      for(PS::S32 j = 0; j < nj; j++) {
	const PS::U32 propj = epj[j].prop;
	const PS::F64vec drij = ri - epj[j].pos;
	const PS::F64 dr2 = drij * drij;
	if(dr2 < Parameter::rc2) {
	  const PS::F64 dr = std::sqrt(dr2);
	  d_sum[propj] += (Parameter::rc - dr) * (Parameter::rc - dr); //NOTE: density kernel is harmonic
	}
      }
      
      for(PS::S32 k = 0; k < Parameter::prop_num; k++)
	result[i].dens[k] += d_sum[k];
    }
  }
};
  
struct CalcForceEpEpDPD {
  //prng seed
  static PS::U32 m_seed;
    
  void operator () (const EPI::DPD* epi,
		    const PS::S32 ni,
		    const EPJ::DPD* epj,
		    const PS::S32 nj,
		    RESULT::ForceDPD* result)
  {
    for(PS::S32 i = 0; i < ni; i++) {
      const PS::F64vec ri = epi[i].pos;
      const PS::F64vec vi = epi[i].vel;
      const PS::U32   idi = epi[i].id;
      const PS::U32 propi = epi[i].prop;
      const std::array<PS::F64, Parameter::prop_num> densi = epi[i].dens;
      
      PS::F64vec fsum(0.0, 0.0, 0.0);
      PS::F64vec psum(0.0, 0.0, 0.0);
      
      for(PS::S32 j = 0; j < nj; j++) {
	if(idi == epj[j].id) continue;
	
	const PS::F64vec drij = ri - epj[j].pos;
	const PS::F64 dr2 = drij * drij;
	if(dr2 < Parameter::rc2) {
	  const PS::F64vec dvij = vi - epj[j].vel;
	  PS::F64 densij[Parameter::prop_num];

	  for(PS::S32 k = 0; k < Parameter::prop_num; k++)
	    densij[k] = densi[k] + epj[j].dens[k];
	  
	  const PS::U32 idj = epj[j].id;
	  
	  PS::U32 m_i = idi, m_j = idj;
	  if(idi > idj){ // m_i <= m_j 
	    m_i = idj;
	    m_j = idi;
	  }

	  Saru saru(m_i, m_j, m_seed + Parameter::time);
	  const PS::F64 rnd = saru.nrml();

	  //kernel for thermostat
	  const PS::F64 dr = std::sqrt(dr2);
	  const PS::F64 inv_dr = 1.0 / dr;
	  const PS::U32 propj  = epj[j].prop;
	  const PS::F64 one_m_dr = 1.0 - dr * Parameter::irc;

	  //kernel for conservative force
#ifdef PAIRWISE_DPD
	  const PS::F64 cf_co = 25.0 * one_m_dr;
	  const PS::F64 cf_mbd = 0.0;
	  assert(std::isfinite(rnd));
#else

	  const PS::F64 cf_co  = Parameter::cf_c[propi][propj] * 6.0 * (dr - Parameter::arc) * (dr - Parameter::rc) * (dr >= Parameter::arc);
	  PS::F64 cf_mbd = 0.0;
	  for(PS::S32 k = 0; k < Parameter::prop_num; k++)
	    cf_mbd += densij[k] * Parameter::cf_m[propi][propj][k];
	  cf_mbd *= one_m_dr;
#endif

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
  
PS::U32 CalcForceEpEpDPD::m_seed = 0;
