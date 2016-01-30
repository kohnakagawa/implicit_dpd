#pragma once

#include "saruprngCUDA.cuh"

__device__ __forceinline__ uint __float_as_uint( float r )
{
  uint u;
  asm volatile( "mov.b32 %0, %1;" : "=r"( u ) : "f"( r ) );
  return u;
}

template<class VecPos, class VecForce, typename T>
__global__ void ForceKernel(const int2* __restrict__ ij_disp,
     	  	            const EPI::DensityGPU<VecPos>* __restrict__ epi,
			    const EPJ::DensityGPU<VecPos>* __restrict__ epj,
			    RESULT::DensityGPU<T>* __restrict__ result,
			    const uint seed) {
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  const VecPos ri = epi[tid].pos;
  T d_sum[Parameter::prop_num] = {(T)0.0};
  
  const int j_head = ij_disp[epi[tid].id_walk    ].y;
  const int j_tail = ij_disp[epi[tid].id_walk + 1].y;
  
  for(int j = j_head; j < j_tail; j++) {
#ifdef USE_FLOAT_VEC
    const VecPos rj = __ldg(&epj[j].pos);
    const uint prpj = __float_as_uint(rj.w);
#else
    const VecPos rj = epj[j].pos;
    const uint prpj = epj[j].prop_;
#endif
    const VecForce drij = {ri.x - rj.x, ri.y - rj.y, ri.z - rj.z};
    const T dr2 = drij.x * drij.x + drij.y * drij.y + drij.z * drij.z;
    const T dr = sqrtf(dr2);
    d_sum[prpj] += (Parameter::rc - dr) * (Parameter::rc - dr) * (dr < Parameter::rc) * (dr != 0.0);
  }

  for(int k = 0; k < Parameter::prop_num; k++)
    result[tid].dens[k] = d_sum[k];
}

template<class VecPos, class VecForce, typename T>
__global__ void ForceKernel(const int2* __restrict__ ij_disp,
			    const EPI::DPDGPU<VecPos>* __restrict__ epi,
			    const EPJ::DPDGPU<VecPos>* __restrict__ epj,
			    RESULT::ForceGPU<VecForce>* __restrict__ force,
			    const uint seed) {
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  
  const VecPos ri = epi[tid].pos;
  const VecPos vi = epi[tid].vel;

#ifdef USE_FLOAT_VEC
  const uint idi   = __float_as_uint(ri.w);
  const uint prpi  = __float_as_uint(vi.w);
#else
  const uint idi   = epi[tid].id_;
  const uint prpi  = epi[tid].prop_;
#endif

  T densi[Parameter::prop_num];
  for(int k = 0; k < Parameter::prop_num; k++)
    densi[k] = epi[tid].dens[k];
  
  const int j_head = ij_disp[epi[tid].id_walk    ].y;
  const int j_tail = ij_disp[epi[tid].id_walk + 1].y;

  VecForce fsum = {(T)0.0};
  VecForce psum = {(T)0.0};
  
  for(int j = j_head; j < j_tail; j++) {
    const VecPos rj = epj[j].pos;
    const VecPos vj = epj[j].vel;

#ifdef USE_FLOAT_VEC    
    const uint idj   = __float_as_uint(rj.w);
    const uint prpj = __float_as_uint(vj.w);
#else
    const uint idj   = epj[j].id_;
    const uint prpj = epj[j].prop_;
#endif    
    
    const VecForce drij = {ri.x - rj.x, ri.y - rj.y, ri.z - rj.z};
    const VecForce dvij = {vi.x - vj.x, vi.y - vj.y, vi.z - vj.z};

    T densij[Parameter::prop_num];
    for(int k = 0; k < Parameter::prop_num; k++)
      densij[k] = densi[k] + epj[j].dens[k];
    
    uint m_i, m_j;
    const bool flag = (idi > idj); //m_i < m_j
    m_i = flag ? idj : idi;
    m_j = flag ? idi : idj;
    
    SaruGPU saru(m_i, m_j, seed);
    const T rnd = saru.nrml_f();
    
    const T dr2 = drij.x * drij.x + drij.y * drij.y + drij.z * drij.z;
    
    //kernel for thermostat
    const T inv_dr = rsqrtf(dr2);
    const T dr = inv_dr * dr2;
    const T one_m_dr = 1.0 - dr * Parameter::irc;

#ifdef PAIRWISE_DPD
    const T cf_co = 25.0 * one_m_dr;
    const T cf_mbd = 0.0;

#else

#ifdef USE_TEXTURE_MEM
    const T cf_co = tex2D(cf_c, prpj, prpi) * (dr - Parameter::arc) * (dr - Parameter::rc) * (dr >= Parameter::arc);

    T cf_mbd = 0.0;
    for(int k = 0; k < Parameter::prop_num; k++)
      cf_mbd += densij[k] * tex3D(cf_m, k, prpj, prpi);
    cf_mbd *= one_m_dr;
    
#else
    const T cf_co = cf_c[prpi][prpj] * (dr - Parameter::arc) * (dr - Parameter::rc) * (dr >= Parameter::arc);

    T cf_mbd = 0.0;
    for(int k = 0; k < Parameter::prop_num; k++)
      cf_mbd += densij[k] * cf_m[prpi][prpj][k];
    cf_mbd *= one_m_dr;
#endif //USE_TEXTURE_MEM

#endif //PAIRWISE_DPD
    const T wrij = one_m_dr;
    //const T wrij = sqrtf(one_m_dr);
    const T sq_wrij = sqrtf(wrij);
    
    const T drij_dvij = drij.x * dvij.x + drij.y * dvij.y + drij.z * dvij.z;

#ifdef USE_TEXTURE_MEM
    T all_cf = (cf_co + cf_mbd +
      		tex2D(cf_r, prpj, prpi) * sq_wrij * rnd -
		tex2D(cf_g, prpj, prpi) * wrij * drij_dvij * inv_dr) * inv_dr;
#else
    T all_cf = (cf_co + cf_mbd +
      		cf_r[prpi][prpj] * sq_wrij * rnd - 
		cf_g[prpi][prpj] * wrij * drij_dvij * inv_dr) * inv_dr;
#endif

    if(dr2 >= Parameter::rc2 || dr2 == 0.0) all_cf = 0.0;
    
    const VecForce dF = {all_cf * drij.x, all_cf * drij.y, all_cf * drij.z};
    const VecForce dP = {dF.x * drij.x, dF.y * drij.y, dF.z * drij.z};
    
    fsum.x += dF.x; fsum.y += dF.y; fsum.z += dF.z;
    psum.x += dP.x; psum.y += dP.y; psum.z += dP.z;
  }
  
  force[tid].acc = fsum;
  force[tid].press = psum;
}
