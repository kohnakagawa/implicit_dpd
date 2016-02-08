#pragma once

#define USE_FLOAT_VEC

#include "f_calculator_gpu.cuh"
#include "ptcl_class.hpp"

bool Policy::Force::init_call = true;
cuda_ptr<EPI::DPDGPU<VecPos> > Policy::Force::dev_epi;
cuda_ptr<EPJ::DPDGPU<VecPos> > Policy::Force::dev_epj;
cuda_ptr<RESULT::ForceGPU<VecForce> > Policy::Force::dev_force;

namespace {
  enum {
    NUM_TOT = 12288,
    //NUM_TOT = 1536,
    PAIR_NUM = 200,
    N_THREAD_GPU = 32,
    N_WALK_LIMIT = 1000,
    NI_LIMIT     = N_WALK_LIMIT*1000,
    NJ_LIMIT     = N_WALK_LIMIT*10000,
  };
  
  bool gpu_inited = false;
  cuda_ptr<int2> ij_disp;
  cuda_ptr<uint> pairs;
  cuda_ptr<uint> num_pairs;
  cuda_ptr<float> length;

}


#include "kernel.cuh"

PS::S32 DispatchTestKernel(const PS::S32 tag,
			   const PS::S32 n_walk,
			   const EPI::DPD ** epi,
			   const PS::S32 * n_epi,
			   const EPJ::DPD ** epj,
			   const PS::S32 * n_epj) {
  assert(n_walk <= N_WALK_LIMIT);

  static int time_cur = 0;

  //allocate array
  if(!gpu_inited) {
    pairs.allocate(NUM_TOT * PAIR_NUM);
    length.allocate(NUM_TOT * PAIR_NUM);
    num_pairs.allocate(NUM_TOT);
    ij_disp.allocate(N_WALK_LIMIT + 2);

    pairs.set_val(0xffffffff);
    length.set_val(0.0);
    num_pairs.set_val(0);

    gpu_inited = true;
  }
  
  if(Policy::Force::init_call) {
    Policy::Force::dev_epi.allocate(NI_LIMIT);
    Policy::Force::dev_epj.allocate(NJ_LIMIT);
    Policy::Force::dev_force.allocate(NI_LIMIT);
    
    Policy::Force::init_call = false;
  }
  
  ij_disp[0].x = ij_disp[0].y = 0;
  for(PS::S32 k = 0; k < n_walk; k++) {
    ij_disp[k + 1].x = ij_disp[k].x + n_epi[k];
    ij_disp[k + 1].y = ij_disp[k].y + n_epj[k];
  }
  ij_disp[n_walk + 1] = ij_disp[n_walk];

  assert(ij_disp[n_walk].x < NI_LIMIT);
  assert(ij_disp[n_walk].y < NJ_LIMIT);
  
  ij_disp.host2dev(0, n_walk + 2);
  
  PS::S32 ni_tot_reg = ij_disp[n_walk].x;
  if(ni_tot_reg % N_THREAD_GPU) {
    ni_tot_reg /= N_THREAD_GPU;
    ni_tot_reg++;
    ni_tot_reg *= N_THREAD_GPU;
  }
  
  PS::S32 ni_tot = -1, nj_tot = -1;
  Policy::Force::CopyToBuffer(n_walk, epi, n_epi, epj, n_epj, ni_tot, nj_tot);
  
  for(PS::S32 i = ni_tot; i < ni_tot_reg; i++)
    Policy::Force::dev_epi[i].id_walk = n_walk;

  Policy::Force::dev_epi.host2dev(0, ni_tot_reg);
  Policy::Force::dev_epj.host2dev(0, nj_tot);
  
  const PS::S32 nblocks = ni_tot_reg / N_THREAD_GPU;
  const PS::S32 nthreads = N_THREAD_GPU;
  
  ForceTestKernel<float4, float4, float> <<< nblocks, nthreads >>> (ij_disp, Policy::Force::dev_epi, Policy::Force::dev_epj, Policy::Force::dev_force, pairs, num_pairs, length, PAIR_NUM, time_cur, ni_tot);
  time_cur++;

  return 0;
}

PS::S32 RetrieveTestKernel(const PS::S32 tag,
			   const PS::S32 n_walk,
			   const PS::S32 ni[],
			   RESULT::ForceDPD * force[]) {
  int ni_tot = 0;
  for(int k = 0; k < n_walk; k++)
    ni_tot += ni[k];
  Policy::Force::dev_force.dev2host(0, ni_tot);
  Policy::Force::CopyToOrigin(n_walk, ni, force);
  return 0;
}
