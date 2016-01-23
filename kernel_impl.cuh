#pragma once

template<class VecPos, class VecForce, typename T>
__global__ void DensityKernel(const int2 *ij_disp,
			      const EPI::DensityGPU<VecPos>* epi,
			      const EPJ::DensityGPU<VecPos>* epj,
			      RESULT::DensityGPU<VecForce>* result) {
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  const VecPos ri = epi[tid].pos;
  T d_sum[Parameter::prop_num] = {0.f};
  
  const int j_head = ij_disp[epi[tid].id_walk    ].y;
  const int j_tail = ij_disp[epi[tid].id_walk + 1].y;
  
  for(int j = j_head; j < j_tail; j++) {
    const VecPos rj = epj[j].pos;
#ifdef USE_FLOAT_VEC
    const uint prpj = __float_as_int(rj.w);
#else
    const uint prpj = epj[j].prop_;
#endif
    const VecForce drij = {ri.x - rj.x, ri.y - rj.y, ri.z - rj.z};
    const T dr2 = drij.x * drij.x + drij.y * drij.y + drij.z * drij.z;
    const T dr = sqrt(dr2);
    d_sum[prpj] += (Parameter::rc - dr) * (Parameter::rc - dr) * (dr < Parameter::rc);
  }

  for(int k = 0; k < Parameter::prop_num; k++)
    result.dens[k] = d_sum[k];
}

template<class VecPos, class VecForce, typename T>
__global__ void ForceKernel(const int2 *ij_disp,
			    const EPI::DPDGPU<VecPos>* epi,
			    const EPJ::DPDGPU<VecPos>* epj,
			    RESULT::ForceGPU<VecForce>* force,
			    const uint seed) {
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  
  const VecPos ri = epi[tid].pos;
  const VecPos vi = epi[tid].vel;

#ifdef USE_FLOAT_VEC
  const uint idi   = __float_as_int(ri.w);
  const uint prpi  = __float_as_int(vi.w);
#else
  const uint idi   = epi[tid].id_;
  const uint prpi  = epi[tid].prop_;
#endif

  T densi[Parameter::prop_num];
  for(int k = 0; k < Parameter::prop_num; k++)
    densi[k] = epi[tid].dens[k];
  
  const int j_head = ij_disp[epi[tid].id_walk    ].y;
  const int j_tail = ij_disp[epi[tid].id_walk + 1].y;
  
  VecForce fsum = {0.f, 0.f, 0.f, 0.f};
  VecForce psum = {0.f, 0.f, 0.f, 0.f};
  for(int j = j_head; j < j_tail; j++) {
    const VecPos rj = epj[j].pos;
    const VecPos vj = epj[j].vel;

#ifdef USE_FLOAT_VEC    
    const uint idj   = __float_as_int(rj.w);
    const uint propj = __float_as_int(vj.w);
#else
    const uint idj   = epj[j].id_;
    const uint propj = epj[j].prop_;
#endif    
    
    const VecForce drij = {ri.x - rj.x, ri.y - rj.y, ri.z - rj.z};
    const VecForce dvij = {vi.x - vj.x, vi.y - vj.y, vi.z - vj.z};

    T densij[Parameter::prop_num];
    for(int k = 0; k < Parameter::prop_num; k++)
      densij[k] = densi[k] + epj[j].dens[k];
    
    uint m_i = idi, m_j = idj;
    if(idi > idj) { //m_i <= m_j
      m_i = idj;
      m_j = idi;
    }
    
    SaruGPU saru(m_i, m_j, seed);
    const T rnd = saru.nrml();
    
    const T dr2 = drij.x * drij.x + drij.y * drij.y + drij.z * drij.z;
    
    //kernel for thermostat
    const T inv_dr = rsqrtf(dr2);
    const T dr = inv_dr * dr2;
    const T one_m_dr = 1.0 - dr * Parameter::irc;

#ifdef PAIRWISE_DPD
    const T cf_co = 25.0 * one_m_dr;
    const T cf_mbd = 0.0;
#else
    const T cf_co = Parameter::cf_c[prpi][prpj] * (dr - Parameter::arc) * (dr - Parameter::rc) * (dr >= Parameter::arc);
    
    const T cf_mbd = 0.0;
    for(int k = 0; k < Parameter::prop_num; k++)
      cf_mbd += densij[k] * const_gpu::cf_m[prpi][prpj][k];
    cf_mbd *= one_m_dr;
#endif
    const T wrij = one_m_dr;
    //const float wrij = sqrtf(one_m_dr);
    const T sq_wrij = sqrtf(wrij);
    
    const T drij_dvij = drij.x * dvij.x + drij.y * dvij.y + drij.z * dvij.z;
    const T all_cf = (cf_co + cf_mbd +
			  const_gpu::cf_r[prpi][prpj] * sq_wrij * rnd - 
			  const_gpu::cf_g[prpi][prpj] * wrij * drij_dvij * inv_dr) * inv_dr;
    if(dr2 < Parameter::dr2) all_cf = 0.0;
    
    const VecForce dF = {all_cf * drij.x, all_cf * drij.y, all_cf * drij.z};
    const VecForce dP = {dF.x * drij.x, dF.y * drij.y, dF.z * drij.z};
    fsum.x += dF.x; fsum.y += dF.y; fsum.z += dF.z;
    psum.x += dP.x; psum.y += dP.y; psum.z += dP.z;
  }
  
  force[tid].acc = fsum;
  force[tid].press = psum;
}
