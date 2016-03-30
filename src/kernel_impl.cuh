#pragma once

#include "saruprngCUDA.cuh"

__device__ __forceinline__ uint __float_as_uint( float r )
{
  uint u;
  asm volatile( "mov.b32 %0, %1;" : "=r"( u ) : "f"( r ) );
  return u;
}

#if 0

template<class VecPos, class VecForce, typename T>
__global__ void ForceKernel(const int2* __restrict__ ij_disp,
     	  	            const EPI::DensityGPU<VecPos>* __restrict__ epi,
			    const EPJ::DensityGPU<VecPos>* __restrict__ epj,
			    RESULT::DensityGPU<T>* __restrict__ result,
			    const uint seed) {
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  const VecPos ri = epi[tid].pos;

  __shared__ T d_sum_sh[N_THREAD_GPU * Parameter::prop_num];
  T* d_sum = &(d_sum_sh[threadIdx.x * Parameter::prop_num]);
#pragma unroll
  for (int k = 0; k < Parameter::prop_num; k++)
    d_sum[k] = 0.0;
  
  const int j_head = ij_disp[epi[tid].id_walk    ].y;
  const int j_tail = ij_disp[epi[tid].id_walk + 1].y;
  
  for (int j = j_head; j < j_tail; j++) {
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
    d_sum[prpj] += ((T)Parameter::rc - dr) * ((T)Parameter::rc - dr) * (dr < (T)Parameter::rc) * (dr != 0.0);
  }

#pragma unroll
  for (int k = 0; k < Parameter::prop_num; k++)
    result[tid].dens[k] = d_sum[k];
}

#else

//tuned version
template<class VecPos, class VecForce, typename T>
__global__ void ForceKernel(const int2* __restrict__ ij_disp,
     	  	            const EPI::DensityGPU<VecPos>* __restrict__ epi,
			    const EPJ::DensityGPU<VecPos>* __restrict__ epj,
			    RESULT::DensityGPU<T>* __restrict__ result,
			    const uint seed) {
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  const VecPos ri = epi[tid].pos;

  __shared__ T d_sum_sh[N_THREAD_GPU * Parameter::prop_num];
  T* d_sum = &(d_sum_sh[threadIdx.x * Parameter::prop_num]);
#pragma unroll
  for (int k = 0; k < Parameter::prop_num; k++)
    d_sum[k] = 0.0;
  
  const int j_head = ij_disp[epi[tid].id_walk    ].y;
  const int j_tail = ij_disp[epi[tid].id_walk + 1].y;
  
  int j = j_head;
  const int res_loop_cnt = ( (j_tail - j_head) & 7);
  const int ini_loop = (res_loop_cnt & 3) + j_head;
  const int nxt_loop = res_loop_cnt + j_head;
  
  for (; j < ini_loop; j++) {
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
    d_sum[prpj] += ((T)Parameter::rc - dr) * ((T)Parameter::rc - dr) * (dr2 < (T)Parameter::rc2) * (dr2 != 0.0);
  } //end of for j

  //unroll 4
  for (; j < nxt_loop; j += 4) {
#ifdef USE_FLOAT_VEC
    const VecPos rj0 = __ldg(&epj[j    ].pos); const uint prpj0 = __float_as_uint(rj0.w);
    const VecPos rj1 = __ldg(&epj[j + 1].pos); const uint prpj1 = __float_as_uint(rj1.w);
    const VecPos rj2 = __ldg(&epj[j + 2].pos); const uint prpj2 = __float_as_uint(rj2.w);
    const VecPos rj3 = __ldg(&epj[j + 3].pos); const uint prpj3 = __float_as_uint(rj3.w);
#else
    const VecPos rj0 = epj[j    ].pos; const uint prpj0 = epj[j    ].prop_;
    const VecPos rj1 = epj[j + 1].pos; const uint prpj1 = epj[j + 1].prop_;
    const VecPos rj2 = epj[j + 2].pos; const uint prpj2 = epj[j + 2].prop_;
    const VecPos rj3 = epj[j + 3].pos; const uint prpj3 = epj[j + 3].prop_;
#endif
    const VecForce drij0 = {ri.x - rj0.x, ri.y - rj0.y, ri.z - rj0.z};
    const VecForce drij1 = {ri.x - rj1.x, ri.y - rj1.y, ri.z - rj1.z};
    const VecForce drij2 = {ri.x - rj2.x, ri.y - rj2.y, ri.z - rj2.z};
    const VecForce drij3 = {ri.x - rj3.x, ri.y - rj3.y, ri.z - rj3.z};

    const T dr2_0 = drij0.x * drij0.x + drij0.y * drij0.y + drij0.z * drij0.z;
    const T dr2_1 = drij1.x * drij1.x + drij1.y * drij1.y + drij1.z * drij1.z;
    const T dr2_2 = drij2.x * drij2.x + drij2.y * drij2.y + drij2.z * drij2.z;
    const T dr2_3 = drij3.x * drij3.x + drij3.y * drij3.y + drij3.z * drij3.z;
    
    const T dr_0 = sqrtf(dr2_0);
    const T dr_1 = sqrtf(dr2_1);
    const T dr_2 = sqrtf(dr2_2);
    const T dr_3 = sqrtf(dr2_3);

    d_sum[prpj0] += ((T)Parameter::rc - dr_0) * ((T)Parameter::rc - dr_0) * (dr2_0 < (T)Parameter::rc2) * (dr2_0 != 0.0);
    d_sum[prpj1] += ((T)Parameter::rc - dr_1) * ((T)Parameter::rc - dr_1) * (dr2_1 < (T)Parameter::rc2) * (dr2_1 != 0.0);
    d_sum[prpj2] += ((T)Parameter::rc - dr_2) * ((T)Parameter::rc - dr_2) * (dr2_2 < (T)Parameter::rc2) * (dr2_2 != 0.0);
    d_sum[prpj3] += ((T)Parameter::rc - dr_3) * ((T)Parameter::rc - dr_3) * (dr2_3 < (T)Parameter::rc2) * (dr2_3 != 0.0);
  }
  
  //unroll 8
  for (; j < j_tail; j += 8) {
#ifdef USE_FLOAT_VEC
    const VecPos rj0 = __ldg(&epj[j    ].pos); const uint prpj0 = __float_as_uint(rj0.w);
    const VecPos rj1 = __ldg(&epj[j + 1].pos); const uint prpj1 = __float_as_uint(rj1.w);
    const VecPos rj2 = __ldg(&epj[j + 2].pos); const uint prpj2 = __float_as_uint(rj2.w);
    const VecPos rj3 = __ldg(&epj[j + 3].pos); const uint prpj3 = __float_as_uint(rj3.w);
    const VecPos rj4 = __ldg(&epj[j + 4].pos); const uint prpj4 = __float_as_uint(rj4.w);
    const VecPos rj5 = __ldg(&epj[j + 5].pos); const uint prpj5 = __float_as_uint(rj5.w);
    const VecPos rj6 = __ldg(&epj[j + 6].pos); const uint prpj6 = __float_as_uint(rj6.w);
    const VecPos rj7 = __ldg(&epj[j + 7].pos); const uint prpj7 = __float_as_uint(rj7.w);
#else
    const VecPos rj0 = epj[j    ].pos; const uint prpj0 = epj[j    ].prop_;
    const VecPos rj1 = epj[j + 1].pos; const uint prpj1 = epj[j + 1].prop_;
    const VecPos rj2 = epj[j + 2].pos; const uint prpj2 = epj[j + 2].prop_;
    const VecPos rj3 = epj[j + 3].pos; const uint prpj3 = epj[j + 3].prop_;
    const VecPos rj4 = epj[j + 4].pos; const uint prpj4 = epj[j + 4].prop_;
    const VecPos rj5 = epj[j + 5].pos; const uint prpj5 = epj[j + 5].prop_;
    const VecPos rj6 = epj[j + 6].pos; const uint prpj6 = epj[j + 6].prop_;
    const VecPos rj7 = epj[j + 7].pos; const uint prpj7 = epj[j + 7].prop_;
#endif
    const VecForce drij0 = {ri.x - rj0.x, ri.y - rj0.y, ri.z - rj0.z};
    const VecForce drij1 = {ri.x - rj1.x, ri.y - rj1.y, ri.z - rj1.z};
    const VecForce drij2 = {ri.x - rj2.x, ri.y - rj2.y, ri.z - rj2.z};
    const VecForce drij3 = {ri.x - rj3.x, ri.y - rj3.y, ri.z - rj3.z};
    const VecForce drij4 = {ri.x - rj4.x, ri.y - rj4.y, ri.z - rj4.z};
    const VecForce drij5 = {ri.x - rj5.x, ri.y - rj5.y, ri.z - rj5.z};
    const VecForce drij6 = {ri.x - rj6.x, ri.y - rj6.y, ri.z - rj6.z};
    const VecForce drij7 = {ri.x - rj7.x, ri.y - rj7.y, ri.z - rj7.z};
    
    const T dr2_0 = drij0.x * drij0.x + drij0.y * drij0.y + drij0.z * drij0.z;
    const T dr2_1 = drij1.x * drij1.x + drij1.y * drij1.y + drij1.z * drij1.z;
    const T dr2_2 = drij2.x * drij2.x + drij2.y * drij2.y + drij2.z * drij2.z;
    const T dr2_3 = drij3.x * drij3.x + drij3.y * drij3.y + drij3.z * drij3.z;
    const T dr2_4 = drij4.x * drij4.x + drij4.y * drij4.y + drij4.z * drij4.z;
    const T dr2_5 = drij5.x * drij5.x + drij5.y * drij5.y + drij5.z * drij5.z;
    const T dr2_6 = drij6.x * drij6.x + drij6.y * drij6.y + drij6.z * drij6.z;
    const T dr2_7 = drij7.x * drij7.x + drij7.y * drij7.y + drij7.z * drij7.z;

    const T dr_0 = sqrtf(dr2_0);
    const T dr_1 = sqrtf(dr2_1);
    const T dr_2 = sqrtf(dr2_2);
    const T dr_3 = sqrtf(dr2_3);
    const T dr_4 = sqrtf(dr2_4);
    const T dr_5 = sqrtf(dr2_5);
    const T dr_6 = sqrtf(dr2_6);
    const T dr_7 = sqrtf(dr2_7);
    
    d_sum[prpj0] += ((T)Parameter::rc - dr_0) * ((T)Parameter::rc - dr_0) * (dr2_0 < (T)Parameter::rc2) * (dr2_0 != 0.0);
    d_sum[prpj1] += ((T)Parameter::rc - dr_1) * ((T)Parameter::rc - dr_1) * (dr2_1 < (T)Parameter::rc2) * (dr2_1 != 0.0);
    d_sum[prpj2] += ((T)Parameter::rc - dr_2) * ((T)Parameter::rc - dr_2) * (dr2_2 < (T)Parameter::rc2) * (dr2_2 != 0.0);
    d_sum[prpj3] += ((T)Parameter::rc - dr_3) * ((T)Parameter::rc - dr_3) * (dr2_3 < (T)Parameter::rc2) * (dr2_3 != 0.0);
    d_sum[prpj4] += ((T)Parameter::rc - dr_4) * ((T)Parameter::rc - dr_4) * (dr2_4 < (T)Parameter::rc2) * (dr2_4 != 0.0);
    d_sum[prpj5] += ((T)Parameter::rc - dr_5) * ((T)Parameter::rc - dr_5) * (dr2_5 < (T)Parameter::rc2) * (dr2_5 != 0.0);
    d_sum[prpj6] += ((T)Parameter::rc - dr_6) * ((T)Parameter::rc - dr_6) * (dr2_6 < (T)Parameter::rc2) * (dr2_6 != 0.0);
    d_sum[prpj7] += ((T)Parameter::rc - dr_7) * ((T)Parameter::rc - dr_7) * (dr2_7 < (T)Parameter::rc2) * (dr2_7 != 0.0);
  } //end of for j

#pragma unroll
  for (int k = 0; k < Parameter::prop_num; k++)
    result[tid].dens[k] = d_sum[k];
}

#endif

#if 0

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
  
  for (int k = 0; k < Parameter::prop_num; k++)
    densi[k] = epi[tid].dens[k];
  
  const int j_head = ij_disp[epi[tid].id_walk    ].y;
  const int j_tail = ij_disp[epi[tid].id_walk + 1].y;

  VecForce fsum = {(T)0.0, (T)0.0, (T)0.0, (T)0.0};
  VecForce psum = {(T)0.0, (T)0.0, (T)0.0, (T)0.0};
  
  for (int j = j_head; j < j_tail; j++) {
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
    for (int k = 0; k < Parameter::prop_num; k++)
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
    const T one_m_dr = 1.0 - dr * (T)Parameter::irc;

#ifdef PAIRWISE_DPD
    const T cf_co = 25.0 * one_m_dr;
    const T cf_mbd = 0.0;
#else
    const T cf_co = tex2D(cf_c, prpj, prpi) * (dr - (T)Parameter::arc) * (dr - (T)Parameter::rc) * (dr >= (T)Parameter::arc);

    T cf_mbd = 0.0;
    for (int k = 0; k < Parameter::prop_num; k++)
      cf_mbd += densij[k] * tex3D(cf_m, k, prpj, prpi);
    cf_mbd *= one_m_dr;
#endif //PAIRWISE_DPD
    
    const T wrij = one_m_dr;
    //const T wrij = sqrtf(one_m_dr);
    const T sq_wrij = sqrtf(wrij);
    
    const T drij_dvij = drij.x * dvij.x + drij.y * dvij.y + drij.z * dvij.z;

    T all_cf = (cf_co + cf_mbd +
      		tex2D(cf_r, prpj, prpi) * sq_wrij * rnd -
		tex2D(cf_g, prpj, prpi) * wrij * drij_dvij * inv_dr) * inv_dr;

    if (dr2 >= (T)Parameter::rc2 || dr2 == 0.0) all_cf = 0.0;
    
    const VecForce dF = {all_cf * drij.x, all_cf * drij.y, all_cf * drij.z};
    const VecForce dP = {dF.x * drij.x, dF.y * drij.y, dF.z * drij.z};
    
    fsum.x += dF.x; fsum.y += dF.y; fsum.z += dF.z;
    psum.x += dP.x; psum.y += dP.y; psum.z += dP.z;
  }
  
  force[tid].acc = fsum;
  force[tid].press = psum;
}

#else

//tuned version
template<class VecPos, class VecForce, typename T>
__global__ void ForceKernel(const int2* __restrict__ ij_disp,
			    const EPI::DPDGPU<VecPos>* __restrict__ epi,
			    const EPJ::DPDGPU<VecPos>* __restrict__ epj,
			    RESULT::ForceGPU<VecForce>* __restrict__ force,
			    const uint seed) {
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  
  const VecPos ri = epi[tid].pos;
  const VecPos vi = epi[tid].vel;

  const uint idi   = __float_as_uint(ri.w);
  const uint prpi  = __float_as_uint(vi.w);

  T densi[Parameter::prop_num]; //NOTE: this is not located at local memory.
  
  for (int k = 0; k < Parameter::prop_num; k++)
    densi[k] = epi[tid].dens[k];
  
  const int j_head = ij_disp[epi[tid].id_walk    ].y;
  const int j_tail = ij_disp[epi[tid].id_walk + 1].y;

  VecForce fsum = {(T)0.0, (T)0.0, (T)0.0, (T)0.0};
  VecForce psum = {(T)0.0, (T)0.0, (T)0.0, (T)0.0};

  int j = j_head;
  const int ini_loop = ((j_tail - j_head) & 3) + j_head;
  for (; j < ini_loop; j++) {
    const VecPos rj = epj[j].pos;
    const VecPos vj = epj[j].vel;

    const uint idj  = __float_as_uint(rj.w);
    const uint prpj = __float_as_uint(vj.w);
    
    const VecForce drij = {ri.x - rj.x, ri.y - rj.y, ri.z - rj.z};
    const VecForce dvij = {vi.x - vj.x, vi.y - vj.y, vi.z - vj.z};

    T densij[Parameter::prop_num];
    for (int k = 0; k < Parameter::prop_num; k++)
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
    const T one_m_dr = 1.0 - dr * (T)Parameter::irc;

    const T cf_co = tex2D(cf_c, prpj, prpi) * (dr - (T)Parameter::arc) * (dr - (T)Parameter::rc) * (dr >= (T)Parameter::arc);

    T cf_mbd = 0.0;
    for (int k = 0; k < Parameter::prop_num; k++)
      cf_mbd += densij[k] * tex3D(cf_m, k, prpj, prpi);
    cf_mbd *= one_m_dr;
    
    const T wrij = one_m_dr;
    //const T wrij = sqrtf(one_m_dr);
    const T sq_wrij = sqrtf(wrij);
    
    const T drij_dvij = drij.x * dvij.x + drij.y * dvij.y + drij.z * dvij.z;

    T all_cf = (cf_co + cf_mbd +
      		tex2D(cf_r, prpj, prpi) * sq_wrij * rnd -
		tex2D(cf_g, prpj, prpi) * wrij * drij_dvij * inv_dr) * inv_dr;

    if(dr2 >= (T)Parameter::rc2 || dr2 == 0.0) all_cf = 0.0;
    
    const VecForce dF = {all_cf * drij.x, all_cf * drij.y, all_cf * drij.z};
    const VecForce dP = {dF.x * drij.x, dF.y * drij.y, dF.z * drij.z};
    
    fsum.x += dF.x; fsum.y += dF.y; fsum.z += dF.z;
    psum.x += dP.x; psum.y += dP.y; psum.z += dP.z;
  } //end of for j
  
  for (; j < j_tail; j += 4) {
    const VecPos rj0 = epj[j    ].pos; const VecPos vj0 = epj[j    ].vel;
    const VecPos rj1 = epj[j + 1].pos; const VecPos vj1 = epj[j + 1].vel;
    const VecPos rj2 = epj[j + 2].pos; const VecPos vj2 = epj[j + 2].vel;
    const VecPos rj3 = epj[j + 3].pos; const VecPos vj3 = epj[j + 3].vel;

    const uint idj0 = __float_as_uint(rj0.w); const uint prpj0 = __float_as_uint(vj0.w);
    const uint idj1 = __float_as_uint(rj1.w); const uint prpj1 = __float_as_uint(vj1.w);
    const uint idj2 = __float_as_uint(rj2.w); const uint prpj2 = __float_as_uint(vj2.w);
    const uint idj3 = __float_as_uint(rj3.w); const uint prpj3 = __float_as_uint(vj3.w);
    
    T densij0[Parameter::prop_num], densij1[Parameter::prop_num], densij2[Parameter::prop_num], densij3[Parameter::prop_num];
#pragma unroll
    for (int k = 0; k < Parameter::prop_num; k++) {
      densij0[k] = densi[k] + epj[j    ].dens[k];
      densij1[k] = densi[k] + epj[j + 1].dens[k];
      densij2[k] = densi[k] + epj[j + 2].dens[k];
      densij3[k] = densi[k] + epj[j + 3].dens[k];
    }
    
    const VecForce drij0 = {ri.x - rj0.x, ri.y - rj0.y, ri.z - rj0.z};
    const VecForce drij1 = {ri.x - rj1.x, ri.y - rj1.y, ri.z - rj1.z};
    const VecForce drij2 = {ri.x - rj2.x, ri.y - rj2.y, ri.z - rj2.z};
    const VecForce drij3 = {ri.x - rj3.x, ri.y - rj3.y, ri.z - rj3.z};

    const VecForce dvij0 = {vi.x - vj0.x, vi.y - vj0.y, vi.z - vj0.z};
    const VecForce dvij1 = {vi.x - vj1.x, vi.y - vj1.y, vi.z - vj1.z};
    const VecForce dvij2 = {vi.x - vj2.x, vi.y - vj2.y, vi.z - vj2.z};
    const VecForce dvij3 = {vi.x - vj3.x, vi.y - vj3.y, vi.z - vj3.z};
    
    uint m_i0, m_j0, m_i1, m_j1, m_i2, m_j2, m_i3, m_j3;
    const bool flag0 = (idi > idj0); //m_i < m_j
    const bool flag1 = (idi > idj1); //m_i < m_j
    const bool flag2 = (idi > idj2); //m_i < m_j
    const bool flag3 = (idi > idj3); //m_i < m_j

    m_i0 = flag0 ? idj0 : idi; m_j0 = flag0 ? idi : idj0;
    m_i1 = flag1 ? idj1 : idi; m_j1 = flag1 ? idi : idj1;
    m_i2 = flag2 ? idj2 : idi; m_j2 = flag2 ? idi : idj2;
    m_i3 = flag3 ? idj3 : idi; m_j3 = flag3 ? idi : idj3;

    SaruGPU saru0(m_i0, m_j0, seed), saru1(m_i1, m_j1, seed), saru2(m_i2, m_j2, seed), saru3(m_i3, m_j3, seed);

    const T dr2_0 = drij0.x * drij0.x + drij0.y * drij0.y + drij0.z * drij0.z;
    const T dr2_1 = drij1.x * drij1.x + drij1.y * drij1.y + drij1.z * drij1.z;
    const T dr2_2 = drij2.x * drij2.x + drij2.y * drij2.y + drij2.z * drij2.z;
    const T dr2_3 = drij3.x * drij3.x + drij3.y * drij3.y + drij3.z * drij3.z;
    
    const T drij_dvij0 = drij0.x * dvij0.x + drij0.y * dvij0.y + drij0.z * dvij0.z;
    const T drij_dvij1 = drij1.x * dvij1.x + drij1.y * dvij1.y + drij1.z * dvij1.z;
    const T drij_dvij2 = drij2.x * dvij2.x + drij2.y * dvij2.y + drij2.z * dvij2.z;
    const T drij_dvij3 = drij3.x * dvij3.x + drij3.y * dvij3.y + drij3.z * dvij3.z;

    const T rnd0 = saru0.nrml_f(), rnd1 = saru1.nrml_f(), rnd2 = saru2.nrml_f(), rnd3 = saru3.nrml_f();
    
    //kernel for thermostat
    const T inv_dr0 = rsqrtf(dr2_0);
    const T inv_dr1 = rsqrtf(dr2_1);
    const T inv_dr2 = rsqrtf(dr2_2);
    const T inv_dr3 = rsqrtf(dr2_3);
    
    const T dr0 = inv_dr0 * dr2_0;
    const T dr1 = inv_dr1 * dr2_1;
    const T dr2 = inv_dr2 * dr2_2;
    const T dr3 = inv_dr3 * dr2_3;
    
    const T one_m_dr0 = 1.0 - dr0 * (T)Parameter::irc;
    const T one_m_dr1 = 1.0 - dr1 * (T)Parameter::irc;
    const T one_m_dr2 = 1.0 - dr2 * (T)Parameter::irc;
    const T one_m_dr3 = 1.0 - dr3 * (T)Parameter::irc;
    
    const T cf_co0 = tex2D(cf_c, prpj0, prpi) * (dr0 - (T)Parameter::arc) * (dr0 - (T)Parameter::rc) * (dr0 >= (T)Parameter::arc);
    const T cf_co1 = tex2D(cf_c, prpj1, prpi) * (dr1 - (T)Parameter::arc) * (dr1 - (T)Parameter::rc) * (dr1 >= (T)Parameter::arc);
    const T cf_co2 = tex2D(cf_c, prpj2, prpi) * (dr2 - (T)Parameter::arc) * (dr2 - (T)Parameter::rc) * (dr2 >= (T)Parameter::arc);
    const T cf_co3 = tex2D(cf_c, prpj3, prpi) * (dr3 - (T)Parameter::arc) * (dr3 - (T)Parameter::rc) * (dr3 >= (T)Parameter::arc);

    T cf_mbd0 = (T)0.0, cf_mbd1 = (T)0.0, cf_mbd2 = (T)0.0, cf_mbd3 = (T)0.0;
#pragma unroll
    for (int k = 0; k < Parameter::prop_num; k++) {
      cf_mbd0 += densij0[k] * tex3D(cf_m, k, prpj0, prpi);
      cf_mbd1 += densij1[k] * tex3D(cf_m, k, prpj1, prpi);
      cf_mbd2 += densij2[k] * tex3D(cf_m, k, prpj2, prpi);
      cf_mbd3 += densij3[k] * tex3D(cf_m, k, prpj3, prpi);
    }
    cf_mbd0 *= one_m_dr0;
    cf_mbd1 *= one_m_dr1;
    cf_mbd2 *= one_m_dr2;
    cf_mbd3 *= one_m_dr3;
    
    const T wrij0 = one_m_dr0;
    const T wrij1 = one_m_dr1;
    const T wrij2 = one_m_dr2;
    const T wrij3 = one_m_dr3;

    // const T wrij0 = sqrtf(one_m_dr0);
    // const T wrij1 = sqrtf(one_m_dr1);
    // const T wrij2 = sqrtf(one_m_dr2);
    // const T wrij3 = sqrtf(one_m_dr3);

    const T sq_wrij0 = sqrtf(wrij0);
    const T sq_wrij1 = sqrtf(wrij1);
    const T sq_wrij2 = sqrtf(wrij2);
    const T sq_wrij3 = sqrtf(wrij3);
    
    T all_cf0 = (cf_co0 + cf_mbd0 +
		 tex2D(cf_r, prpj0, prpi) * sq_wrij0 * rnd0 -
		 tex2D(cf_g, prpj0, prpi) * wrij0 * drij_dvij0 * inv_dr0) * inv_dr0;
    T all_cf1 = (cf_co1 + cf_mbd1 +
		 tex2D(cf_r, prpj1, prpi) * sq_wrij1 * rnd1 -
		 tex2D(cf_g, prpj1, prpi) * wrij1 * drij_dvij1 * inv_dr1) * inv_dr1;
    T all_cf2 = (cf_co2 + cf_mbd2 +
		 tex2D(cf_r, prpj2, prpi) * sq_wrij2 * rnd2 -
		 tex2D(cf_g, prpj2, prpi) * wrij2 * drij_dvij2 * inv_dr2) * inv_dr2;
    T all_cf3 = (cf_co3 + cf_mbd3 +
		 tex2D(cf_r, prpj3, prpi) * sq_wrij3 * rnd3 -
		 tex2D(cf_g, prpj3, prpi) * wrij3 * drij_dvij3 * inv_dr3) * inv_dr3;

    if(dr2_0 >= (T)Parameter::rc2 || dr2_0 == 0.0) all_cf0 = 0.0;
    if(dr2_1 >= (T)Parameter::rc2 || dr2_1 == 0.0) all_cf1 = 0.0;
    if(dr2_2 >= (T)Parameter::rc2 || dr2_2 == 0.0) all_cf2 = 0.0;
    if(dr2_3 >= (T)Parameter::rc2 || dr2_3 == 0.0) all_cf3 = 0.0;

    const VecForce dF0 = {all_cf0 * drij0.x, all_cf0 * drij0.y, all_cf0 * drij0.z};
    const VecForce dF1 = {all_cf1 * drij1.x, all_cf1 * drij1.y, all_cf1 * drij1.z};
    const VecForce dF2 = {all_cf2 * drij2.x, all_cf2 * drij2.y, all_cf2 * drij2.z};
    const VecForce dF3 = {all_cf3 * drij3.x, all_cf3 * drij3.y, all_cf3 * drij3.z};
    
    const VecForce dP0 = {dF0.x * drij0.x, dF0.y * drij0.y, dF0.z * drij0.z};
    const VecForce dP1 = {dF1.x * drij1.x, dF1.y * drij1.y, dF1.z * drij1.z};
    const VecForce dP2 = {dF2.x * drij2.x, dF2.y * drij2.y, dF2.z * drij2.z};
    const VecForce dP3 = {dF3.x * drij3.x, dF3.y * drij3.y, dF3.z * drij3.z};

    fsum.x += dF0.x; fsum.y += dF0.y; fsum.z += dF0.z;
    fsum.x += dF1.x; fsum.y += dF1.y; fsum.z += dF1.z;
    fsum.x += dF2.x; fsum.y += dF2.y; fsum.z += dF2.z;
    fsum.x += dF3.x; fsum.y += dF3.y; fsum.z += dF3.z;

    psum.x += dP0.x; psum.y += dP0.y; psum.z += dP0.z;
    psum.x += dP1.x; psum.y += dP1.y; psum.z += dP1.z;
    psum.x += dP2.x; psum.y += dP2.y; psum.z += dP2.z;
    psum.x += dP3.x; psum.y += dP3.y; psum.z += dP3.z;
  } //end of for j
  
  force[tid].acc = fsum;
  force[tid].press = psum;
}

#endif