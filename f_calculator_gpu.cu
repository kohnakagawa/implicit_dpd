#include "f_calculator_gpu.cuh"
#include "ptcl_class.hpp"

PS::F64 Parameter::cf_c[Parameter::prop_num][Parameter::prop_num];
PS::F64 Parameter::cf_g[Parameter::prop_num][Parameter::prop_num];
PS::F64 Parameter::cf_r[Parameter::prop_num][Parameter::prop_num];
PS::F64 Parameter::cf_m[Parameter::prop_num][Parameter::prop_num][Parameter::prop_num];

bool DensityPolicy::init_call = true;
cuda_ptr<EPI::DensityGPU<VecPos> > DensityPolicy::dev_epi;
cuda_ptr<EPJ::DensityGPU<VecPos> > DensityPolicy::dev_epj;
cuda_ptr<RESULT::DensityGPU<Dtype> > DensityPolicy::dev_force;

bool ForcePolicy::init_call = true;
cuda_ptr<EPI::DPDGPU<VecPos> > ForcePolicy::dev_epi;
cuda_ptr<EPJ::DPDGPU<VecPos> > ForcePolicy::dev_epj;
cuda_ptr<RESULT::ForceGPU<VecForce> > ForcePolicy::dev_force;

namespace {
  __device__ float cf_m[Parameter::prop_num][Parameter::prop_num][Parameter::prop_num];
  __device__ float cf_c[Parameter::prop_num][Parameter::prop_num];
  __device__ float cf_r[Parameter::prop_num][Parameter::prop_num];
  __device__ float cf_g[Parameter::prop_num][Parameter::prop_num];

  enum {
    N_THREAD_GPU = 32,
    N_WALK_LIMIT = 1000,
    NI_LIMIT     = N_WALK_LIMIT*1000,
    NJ_LIMIT     = N_WALK_LIMIT*10000,
  };

  bool gpu_inited = false;
  cuda_ptr<int2> ij_disp;

  void set_const_gpu() {
    float cf_m_h[Parameter::prop_num][Parameter::prop_num][Parameter::prop_num];
    float cf_c_h[Parameter::prop_num][Parameter::prop_num];
    float cf_r_h[Parameter::prop_num][Parameter::prop_num];
    float cf_g_h[Parameter::prop_num][Parameter::prop_num];
    
    for(int i = 0; i < Parameter::prop_num; i++) {
      for(int j = 0; j < Parameter::prop_num; j++) {
	cf_c_h[i][j] = Parameter::cf_c[i][j];
	cf_r_h[i][j] = Parameter::cf_r[i][j];
	cf_g_h[i][j] = Parameter::cf_g[i][j];
      }
    }

    for(int i = 0; i < Parameter::prop_num; i++)
      for(int j = 0; j < Parameter::prop_num; j++)
	for(int k = 0; k < Parameter::prop_num; k++)
	  cf_m_h[i][j][k] = Parameter::cf_m[i][j][k];
    
    const size_t cf_m_size = sizeof(float) * Parameter::prop_num * Parameter::prop_num * Parameter::prop_num;
    const size_t cf_p_size = sizeof(float) * Parameter::prop_num * Parameter::prop_num;
    void* ptr_dev = NULL;

#define COPY_TO_DEVICE_SYMBOL(dev, host, size) \
    checkCudaErrors(cudaGetSymbolAddress(&ptr_dev, dev)); \
    checkCudaErrors(cudaMemcpy(ptr_dev, host, size, cudaMemcpyHostToDevice))
    
    COPY_TO_DEVICE_SYMBOL(cf_m, cf_m_h, cf_m_size);
    COPY_TO_DEVICE_SYMBOL(cf_c, cf_c_h, cf_p_size);
    COPY_TO_DEVICE_SYMBOL(cf_r, cf_r_h, cf_p_size);
    COPY_TO_DEVICE_SYMBOL(cf_g, cf_g_h, cf_p_size);

#undef COPY_TO_DEVICE_SYMBOL
  }
};

#include "kernel_impl.cuh"

template<class Policy, class EPI, class EPJ>
PS::S32 DispatchKernel(const PS::S32 tag,
		       const PS::S32 n_walk,
		       const EPI ** epi,
		       const PS::S32 * n_epi,
		       const EPJ ** epj,
		       const PS::S32 * n_epj) {
  assert(n_walk <= N_WALK_LIMIT);

  //allocate array
  if(!gpu_inited) {
    set_const_gpu();
    ij_disp.allocate(N_WALK_LIMIT + 2);    
    gpu_inited = true;
  }
  
  if(Policy::init_call) {
    Policy::dev_epi.allocate(NI_LIMIT);
    Policy::dev_epj.allocate(NJ_LIMIT);
    Policy::dev_force.allocate(NI_LIMIT);
    
    Policy::init_call = false;
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
  Policy::CopyToBuffer(n_walk, epi, n_epi, epj, n_epj, ni_tot, nj_tot);
  
  for(PS::S32 i = ni_tot; i < ni_tot_reg; i++)
    Policy::dev_epi[i].id_walk = n_walk;

  Policy::dev_epi.host2dev(0, ni_tot_reg);
  Policy::dev_epj.host2dev(0, nj_tot);
  
  const PS::S32 nblocks = ni_tot_reg / N_THREAD_GPU;
  const PS::S32 nthreads = N_THREAD_GPU;
  ForceKernel<VecPos, VecForce, Dtype> <<< nblocks, nthreads >>> (ij_disp, Policy::dev_epi, Policy::dev_epj, Policy::dev_force, Parameter::time);

  return 0;
}

template<class Policy, class RESULT>
PS::S32 RetrieveKernel(const PS::S32 tag,
		       const PS::S32 n_walk,
		       const PS::S32 ni[],
		       RESULT * force[]) {
  int ni_tot = 0;
  for(int k = 0; k < n_walk; k++)
    ni_tot += ni[k];
  Policy::dev_force.dev2host(0, ni_tot);
  Policy::CopyToOrigin(n_walk, ni, force);
  return 0;
}

//instantiate
template
PS::S32 DispatchKernel<DensityPolicy, EPI::Density, EPJ::Density>(const PS::S32 tag,
								  const PS::S32 n_walk,
								  const EPI::Density ** epi,
								  const PS::S32 * n_epi,
								  const EPJ::Density ** epj,
								  const PS::S32 * n_epj);

template
PS::S32 DispatchKernel<ForcePolicy,   EPI::DPD,     EPJ::DPD    >(const PS::S32 tag,
								  const PS::S32 n_walk,
								  const EPI::DPD ** epi,
								  const PS::S32 * n_epi,
								  const EPJ::DPD ** epj,
								  const PS::S32 * n_epj);

template
PS::S32 RetrieveKernel<DensityPolicy, RESULT::Density> (const PS::S32 tag,
					 const PS::S32 n_walk,
					 const PS::S32 ni[],
					 RESULT::Density * force[]);

template
PS::S32 RetrieveKernel<ForcePolicy, RESULT::ForceDPD> (const PS::S32 tag,
					  const PS::S32 n_walk,
					  const PS::S32 ni[],
					  RESULT::ForceDPD * force[]);
