#pragma once

#include "saruprng.hpp"
#include "ptcl_class.hpp"

PS::F64 Parameter::cf_c[Parameter::prop_num][Parameter::prop_num];
PS::F64 Parameter::cf_g[Parameter::prop_num][Parameter::prop_num];
PS::F64 Parameter::cf_r[Parameter::prop_num][Parameter::prop_num];
PS::F64 Parameter::cf_m[Parameter::prop_num][Parameter::prop_num][Parameter::prop_num];

struct CalcDensity {
  void operator () (const EPI::Density* __restrict epi,
                    const PS::S32 ni,
                    const EPJ::Density* __restrict epj,
                    const PS::S32 nj,
                    RESULT::Density* __restrict result)
  {
    for (PS::S32 i = 0; i < ni; i++) {
      const PS::F64vec ri = epi[i].pos;
      PS::F64 d_sum[Parameter::prop_num] {};

      // NOTE: cm_pos is only used when i particle's property is hydrophlic.
      PS::F64vec hypb_pos_sum(0.0, 0.0, 0.0), hypl_pos_sum = ri;
      PS::U32 hypb_cnt = 0, hypl_cnt = 1;

      for (PS::S32 j = 0; j < nj; j++) {
        const PS::U32 propj   = epj[j].prop;
        const PS::F64vec drij = ri - epj[j].pos;
        const PS::F64 dr2     = drij * drij;

        if (dr2 < Parameter::rn_c2) {
          // calc cmpos of hyphob
          if (propj == Parameter::Hyphob) {
            hypb_pos_sum += epj[j].pos;
            hypb_cnt++;
          }

          if (dr2 < 1.0) {
            // density calc
            const PS::F64 dr = std::sqrt(dr2);
            d_sum[propj] += (1.0 - dr) * (1.0 - dr);

            // calc cmpos of hyphil
            if (propj == Parameter::Hyphil) {
              hypl_pos_sum += epj[j].pos;
              hypl_cnt++;
            }
          }
        }
      }

      for (PS::S32 k = 0; k < Parameter::prop_num; k++) {
        result[i].dens[k] += d_sum[k];
      }
      result[i].nei_pos_sum[Parameter::Hyphil] += hypl_pos_sum;
      result[i].nei_cnt[Parameter::Hyphil] += hypl_cnt;
      result[i].nei_pos_sum[Parameter::Hyphob] += hypb_pos_sum;
      result[i].nei_cnt[Parameter::Hyphob] += hypb_cnt;
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
    for (PS::S32 i = 0; i < ni; i++) {
      const PS::F64vec ri = epi[i].pos;
      const PS::F64vec vi = epi[i].vel;
      const PS::U32   idi = epi[i].id;
      const PS::U32 propi = epi[i].prop;
      const std::array<PS::F64, Parameter::prop_num> densi = epi[i].dens;

      PS::F64vec fsum(0.0, 0.0, 0.0);
      PS::F64vec psum(0.0, 0.0, 0.0);

      for (PS::S32 j = 0; j < nj; j++) {
        if (idi == epj[j].id) continue;

        const PS::F64vec drij = ri - epj[j].pos;
        const PS::F64    dr2  = drij * drij;
        if (dr2 < 1.0) {
          const PS::F64vec dvij = vi - epj[j].vel;
          const PS::U32    idj  = epj[j].id;

          // density term
          PS::F64 densij[Parameter::prop_num];
          for (PS::S32 k = 0; k < Parameter::prop_num; k++) {
            densij[k] = densi[k] + epj[j].dens[k];
          }

          // get random noise
          PS::U32 m_i = idi, m_j = idj;
          if (idi > idj){
            m_i = idj;
            m_j = idi;
          }
          Saru saru(m_i, m_j, m_seed + Parameter::time);
          const PS::F64 rnd = saru.nrml();

          //kernel for thermostat
          const PS::F64 dr      = std::sqrt(dr2);
          const PS::F64 inv_dr  = 1.0 / dr;
          const PS::U32 propj   = epj[j].prop;
          const PS::F64 rc_m_dr = 1.0 - dr;

          //kernel for conservative force
#ifdef PAIRWISE_DPD
          (void) densij;
          const PS::F64 cf_pair = 25.0;
          const PS::F64 cf_mbd  = 0.0;
#else
          const PS::F64 cf_pair = Parameter::cf_c[propi][propj] * (dr - Parameter::arc) * (dr >= Parameter::arc);
          PS::F64 cf_mbd = 0.0;
          for (PS::S32 k = 0; k < Parameter::prop_num; k++) {
            cf_mbd += Parameter::cf_m[propi][propj][k] * densij[k];
          }
#endif
          const PS::F64 wrij    = rc_m_dr; // pow = 1
          //const PS::F64 wrij = std::sqrt(rc_m_dr); //pow = 1 / 2
          const PS::F64 sq_wrij = std::sqrt(wrij);

          const PS::F64 drij_dvij = drij * dvij;
          const PS::F64 all_cf = ((cf_pair + cf_mbd) * rc_m_dr +  //conservative
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
