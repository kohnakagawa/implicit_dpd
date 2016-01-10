#pragma once

#include "saruprng.hpp"
// #include "static_for.hpp"
#include "epiepj_tune.hpp"

struct NearlistBuffer {
  PS::S32 i_buf_size = -1;
  static PS::ReallocatableArray<PS::S32>* neigh_id; //[pid][pjd]
  enum {
    J_BUFFER_SIZE = 500,
  };
  
  NearlistBuffer() {}
  ~NearlistBuffer() {
    delete [] neigh_id;
  }
  
  void Initialize(const PS::S32 i_buf_size_) {
    i_buf_size = i_buf_size_;
    neigh_id = new PS::ReallocatableArray<PS::S32> [i_buf_size];
    for(PS::S32 i = 0; i < i_buf_size; i++)
      neigh_id[i].reserve(J_BUFFER_SIZE);
  }
  
  void ClearList() {
    for(PS::S32 i = 0; i < i_buf_size; ++i)
      neigh_id[i].clearSize();
  }
  
  void CheckBufferSize(const PS::S32 num_loc) {
    if(num_loc >= i_buf_size) {
      std::cerr << "increase buffer size\n";
      delete [] neigh_id;
      const PS::S32 new_size = num_loc * 10;
      neigh_id = new PS::ReallocatableArray<PS::S32> [new_size];
      i_buf_size = new_size;
    }
  }
  
};
PS::ReallocatableArray<PS::S32>* NearlistBuffer::neigh_id = nullptr;

struct CalcDensity {
  void operator () (const EPI::DPD_Density* __restrict epi,
		    const PS::S32 ni,
		    const EPJ::DPD_Density* __restrict epj,
		    const PS::S32 nj,
		    RESULT::ForceDPDDensity* __restrict result)
  {
    for(PS::S32 i = 0; i < ni; i++) {
      const PS::F64vec ri = epi[i].pos;
      const PS::U32 idi_loc = epi[i].id_loc;
      PS::F64 d_sum[Parameter::prop_num] = { 0.0 };
      
      for(PS::S32 j = 0; j < nj; j++) {
	const PS::U32 propj = epj[j].prop;
	const PS::F64vec drij = ri - epj[j].pos;
	const PS::F64 dr2 = drij * drij;
	if(dr2 < Parameter::rc2 && dr2 != 0.0) {
	  const PS::F64 dr = std::sqrt(dr2);
	  d_sum[propj] += (Parameter::rc - dr) * (Parameter::rc - dr);
	  NearlistBuffer::neigh_id[idi_loc].push_back(j); //NOTE: push_back may be replaced by pushBackNoCheck
	}
      }
      
      for(PS::S32 k = 0; k < Parameter::prop_num; k++)
	result[i].dens[k] += d_sum[k];
    }
  }
};

/* old version
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
      
      PS::F64vec drija = ri - epj[0].pos;
      PS::U32 prpja = epj[0].prop, prpjb = 0;
      PS::F64 dw = 0.0;
      
      for(PS::S32 j = 1; j < nj + 1; j++) {
	const PS::F64vec drij = drija;
	const PS::F64 dr2 = drij * drij;
	const PS::U32 prpj = prpja;
	prpja = epj[j].prop; //NOTE: In case of j == nj, invalid read occurs.
	drija = ri - epj[j].pos;

	if(dr2 >= Parameter::rc2 || dr2 == 0.0) continue;
	
	//store previous result
	d_sum[prpjb] += dw;
	
	//calculate next weight
	const PS::F64 dr = std::sqrt(dr2);
	dw = (Parameter::rc - dr) * (Parameter::rc - dr);
	prpjb = prpj;
      }
      d_sum[prpjb] += dw;
      
      for(PS::S32 k = 0; k < Parameter::prop_num; k++)
	result[i].dens[k] += d_sum[k];
    }
  }
  };*/

struct CalcForceEpEpDPD {
  //prng seed
  static PS::U32 m_seed;
    
  void operator () (const EPI::DPD_Density* __restrict epi,
		    const PS::S32 ni,
		    const EPJ::DPD_Density* __restrict epj,
		    const PS::S32 nj,
		    RESULT::ForceDPDDensity* __restrict result)
  {
    for(PS::S32 i = 0; i < ni; i++) {
      const PS::F64vec ri	= epi[i].pos;
      const PS::F64vec vi	= epi[i].vel;
      const PS::U32   idi	= epi[i].id;
      const PS::U32 idi_loc	= epi[i].id_loc;
      const PS::U32 propi	= epi[i].prop;
      const std::array<PS::F64, Parameter::prop_num> densi = epi[i].dens;
      
      PS::F64vec fsum(0.0, 0.0, 0.0);
      PS::F64vec psum(0.0, 0.0, 0.0);

      const PS::S32 num_intr = NearlistBuffer::neigh_id[idi_loc].size();
      for(PS::S32 jj = 0; jj < num_intr; jj++) {
	const PS::S32 j = NearlistBuffer::neigh_id[idi_loc][jj];
	const PS::F64vec drij = ri - epj[j].pos;
	const PS::F64 dr2 = drij * drij;
	const PS::F64 dr = std::sqrt(dr2);

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
      result[i].acc += fsum;
      result[i].press += psum;
    }
  }
};
  
PS::U32 CalcForceEpEpDPD::m_seed = 0;
