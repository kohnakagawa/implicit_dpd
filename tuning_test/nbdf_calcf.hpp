#pragma once

#include "saruprng.hpp"
// #include "static_for.hpp"
#include "epiepj.hpp"

struct CalcDensity {
  void operator () (const EPI::Density* __restrict epi,
		    const PS::S32 ni,
		    const EPJ::Density* __restrict epj,
		    const PS::S32 nj,
		    RESULT::Density* __restrict result)
  {
    for(PS::S32 i = 0; i < ni; i++) {
      const PS::F64vec ri = epi[i].pos;
      PS::F64 d_sum[Parameter::prop_num] = { 0.0 };
      for(PS::S32 j = 0; j < nj; j++) {
	const PS::U32 propj = epj[j].prop;
	const PS::F64vec drij = ri - epj[j].pos;
	const PS::F64 dr2 = drij * drij;
	if(dr2 < Parameter::rc2 && dr2 != 0.0) {
	  const PS::F64 dr = std::sqrt(dr2);
	  d_sum[propj] += (Parameter::rc - dr) * (Parameter::rc - dr);
	}
      }
      
      for(PS::S32 k = 0; k < Parameter::prop_num; k++)
	result[i].dens[k] += d_sum[k];
    }
  }
};

struct CalcDensityUnroll {
  void operator () (const EPI::Density* __restrict epi,
		    const PS::S32 ni,
		    const EPJ::Density* __restrict epj,
		    const PS::S32 nj,
		    RESULT::Density* __restrict result)
  {
    for(PS::S32 i = 0; i < ni; i++) {
      const PS::F64vec ri = epi[i].pos;
      PS::F64 d_sum[Parameter::prop_num] = { 0.0 };
      
      PS::S32 j;
      for(j = 0; j < nj; j += 4) {
	const PS::U32 propja = epj[j    ].prop;
	const PS::U32 propjb = epj[j + 1].prop;
	const PS::U32 propjc = epj[j + 2].prop;
	const PS::U32 propjd = epj[j + 3].prop;
	
	const PS::F64vec drija = ri - epj[j    ].pos;
	const PS::F64vec drijb = ri - epj[j + 1].pos;
	const PS::F64vec drijc = ri - epj[j + 2].pos;
	const PS::F64vec drijd = ri - epj[j + 3].pos;
	
	const PS::F64 dr2a = drija * drija;
	const PS::F64 dr2b = drijb * drijb;
	const PS::F64 dr2c = drijc * drijc;
	const PS::F64 dr2d = drijd * drijd;
	
	if(dr2a < Parameter::rc2 && dr2a != 0.0) {
	  const PS::F64 dr = std::sqrt(dr2a);
	  d_sum[propja] += (Parameter::rc - dr) * (Parameter::rc - dr);
	}

	if(dr2b < Parameter::rc2 && dr2b != 0.0) {
	  const PS::F64 dr = std::sqrt(dr2b);
	  d_sum[propjb] += (Parameter::rc - dr) * (Parameter::rc - dr);
	}

	if(dr2c < Parameter::rc2 && dr2c != 0.0) {
	  const PS::F64 dr = std::sqrt(dr2c);
	  d_sum[propjc] += (Parameter::rc - dr) * (Parameter::rc - dr);
	}

	if(dr2d < Parameter::rc2 && dr2d != 0.0) {
	  const PS::F64 dr = std::sqrt(dr2d);
	  d_sum[propjd] += (Parameter::rc - dr) * (Parameter::rc - dr);
	}
      }
      
      j -= 4;
      for(; j < nj; j++) {
	const PS::U32 propj = epj[j].prop;
	const PS::F64vec drij = ri - epj[j].pos;
	const PS::F64 dr2 = drij * drij;
	if(dr2 < Parameter::rc2 && dr2 != 0.0) {
	  const PS::F64 dr = std::sqrt(dr2);
	  d_sum[propj] += (Parameter::rc - dr) * (Parameter::rc - dr);
	}
      }
      
      for(PS::S32 k = 0; k < Parameter::prop_num; k++)
	result[i].dens[k] += d_sum[k];
    }
  }
};

struct CalcDensitySoftPipe {
  void operator () (const EPI::Density* __restrict epi,
		    const PS::S32 ni,
		    const EPJ::Density* __restrict epj,
		    const PS::S32 nj,
		    RESULT::Density* __restrict result)
  {
    for(PS::S32 i = 0; i < ni; i++) {
      const PS::F64vec ri = epi[i].pos;
      PS::F64 d_sum[Parameter::prop_num] = { 0.0 };
      
      PS::F64vec drija = ri - epj[0].pos, drijb;
      PS::F64 dr2a = 0.0;
      PS::U32 prpja = epj[0].prop, prpjb = 0;
      
      for(PS::S32 j = 1; j < nj; j++) {
	dr2a = drija * drija;
	drijb = ri - epj[j].pos;
	prpjb = epj[j].prop;
	
	if(dr2a < Parameter::rc2 && dr2a != 0.0) {
	  const PS::F64 dr = std::sqrt(dr2a);
	  d_sum[prpja] += (Parameter::rc - dr) * (Parameter::rc - dr);
	}
	drija = drijb;
	prpja = prpjb;
      }
      dr2a = drija * drija;
      if(dr2a < Parameter::rc2 && dr2a != 0.0) {
	const PS::F64 dr = std::sqrt(dr2a);
	d_sum[prpja] += (Parameter::rc - dr) * (Parameter::rc - dr);
      }
      
      for(PS::S32 k = 0; k < Parameter::prop_num; k++)
	result[i].dens[k] += d_sum[k];
    }
  }
};

struct CalcForceEpEpDPD {
  //prng seed
  static PS::U32 m_seed;
    
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

#pragma novector
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

	  const PS::F64 cf_co  = Parameter::cf_c[propi][propj] * (dr - Parameter::arc) * (dr - Parameter::rc) * (dr >= Parameter::arc);
	  PS::F64 cf_mbd = 0.0;
#pragma novector
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
