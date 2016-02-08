#pragma once

__device__ __forceinline__ uint __float_as_uint( float r )
{
  uint u;
  asm volatile( "mov.b32 %0, %1;" : "=r"( u ) : "f"( r ) );
  return u;
}

template<class VecPos, class VecForce, typename T>
__global__ void ForceTestKernel(const int2*  ij_disp,
				const EPI::DPDGPU<VecPos>*  epi,
				const EPJ::DPDGPU<VecPos>*  epj,
				RESULT::ForceGPU<VecForce>*  force,
				uint*  pairs,
				uint*  num_pairs,
				float*  length,
				const int buf_size,
				const int time,
				const int real_size) {
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  
  if (tid < real_size) {
    const VecPos ri = epi[tid].pos;
#ifdef USE_FLOAT_VEC
    const uint idi   = __float_as_uint(ri.w);
#else
    const uint idi   = epi[tid].id_;
#endif

    const int j_head = ij_disp[epi[tid].id_walk    ].y;
    const int j_tail = ij_disp[epi[tid].id_walk + 1].y;

    uint* beg_pairs = &(pairs[idi * buf_size]);
    float* beg_leng = &(length[idi * buf_size]);

    uint cnt = num_pairs[idi];
    if(cnt == 0) {
      beg_pairs[0] = idi;
      cnt++;
    }
  
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
      const T dr2 = drij.x * drij.x + drij.y * drij.y + drij.z * drij.z;

      if((dr2 < rc2) && (idi != idj) && (tid < real_size)) {
	beg_pairs[cnt] = idj;
	beg_leng[cnt]  = dr2;
	cnt++;
      }
    }
    num_pairs[idi] = cnt;
  }
}
