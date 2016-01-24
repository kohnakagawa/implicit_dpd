#pragma once

#include "particle_simulator.hpp"
#include "cuda_ptr.cuh"
#include "user_defs.h"
#include "ptcl_class_gpu.cuh"
#include "parameter.hpp"

template<class Policy, class EPI, class EPJ>
PS::S32 DispatchKernel(const PS::S32 tag,
		       const PS::S32 n_walk,
		       const EPI ** epi,
		       const PS::S32 * n_epi,
		       const EPJ ** epj,
		       const PS::S32 * n_epj);

template<class Policy, class RESULT>
PS::S32 RetrieveKernel(const PS::S32 tag,
	     	       const PS::S32 n_walk,
		       const PS::S32 ni[],
		       RESULT * force[]);

void cleanup_GPU();

struct DensityPolicy {
  static bool init_call;
  static cuda_ptr<EPI::DensityGPU<VecPos> > dev_epi;
  static cuda_ptr<EPJ::DensityGPU<VecPos> > dev_epj;
  static cuda_ptr<RESULT::DensityGPU<Dtype> > dev_force;
  
  template<class EPI, class EPJ>
  static void CopyToBuffer(const PS::S32 n_walk,
			   const EPI** epi,
			   const PS::S32* n_epi,
			   const EPJ** epj,
			   const PS::S32* n_epj,
			   PS::S32& ni_tot,
			   PS::S32& nj_tot) {
    ni_tot = nj_tot = 0;
    for(PS::S32 iw = 0; iw < n_walk; iw++) {
      for(PS::S32 i = 0; i < n_epi[iw]; i++) {
	dev_epi[ni_tot].pos.x = epi[iw][i].pos.x;
	dev_epi[ni_tot].pos.y = epi[iw][i].pos.y;
	dev_epi[ni_tot].pos.z = epi[iw][i].pos.z;

	dev_epi[ni_tot].prop() = epi[iw][i].prop;

	dev_epi[ni_tot].id_walk = iw;
	ni_tot++;
      }
    
      for(PS::S32 j = 0; j < n_epj[iw]; j++) {
	dev_epj[nj_tot].pos.x = epj[iw][j].pos.x;
	dev_epj[nj_tot].pos.y = epj[iw][j].pos.y;
	dev_epj[nj_tot].pos.z = epj[iw][j].pos.z;

	dev_epj[nj_tot].prop() = epj[iw][j].prop;
      
	nj_tot++;
      }
    }
  }
  
  template<class RESULT>
  static void CopyToOrigin(const PS::S32 n_walk,
			   const PS::S32 ni[],
			   RESULT* force[]) {
    int n_cnt = 0;
    for(int iw = 0; iw < n_walk; iw++) {
      for(int i = 0; i < ni[iw]; i++) {

	for(int k = 0; k < Parameter::prop_num; k++)
	  force[iw][i].dens[k] = dev_force[n_cnt].dens[k];

	n_cnt++;
      }
    }
  }
};

struct ForcePolicy {
  static bool init_call;
  static cuda_ptr<EPI::DPDGPU<VecPos> > dev_epi;
  static cuda_ptr<EPJ::DPDGPU<VecPos> > dev_epj;
  static cuda_ptr<RESULT::ForceGPU<VecForce> > dev_force;
  
  template<class EPI, class EPJ>
  static void CopyToBuffer(const PS::S32 n_walk,
			   const EPI** epi,
			   const PS::S32* n_epi,
			   const EPJ** epj,
			   const PS::S32* n_epj,
			   PS::S32& ni_tot,
			   PS::S32& nj_tot) {
    ni_tot = nj_tot = 0;
    for(PS::S32 iw = 0; iw < n_walk; iw++) {
      for(PS::S32 i = 0; i < n_epi[iw]; i++) {
	dev_epi[ni_tot].pos.x = epi[iw][i].pos.x;
	dev_epi[ni_tot].pos.y = epi[iw][i].pos.y;
	dev_epi[ni_tot].pos.z = epi[iw][i].pos.z;
      
	dev_epi[ni_tot].vel.x = epi[iw][i].vel.x;
	dev_epi[ni_tot].vel.y = epi[iw][i].vel.y;
	dev_epi[ni_tot].vel.z = epi[iw][i].vel.z;
      
	dev_epi[ni_tot].id()  = epi[iw][i].id;
	dev_epi[ni_tot].prop() = epi[iw][i].prop;

	dev_epi[ni_tot].id_walk = iw;
	ni_tot++;
      }
    
      for(PS::S32 j = 0; j < n_epj[iw]; j++) {
	dev_epj[nj_tot].pos.x = epj[iw][j].pos.x;
	dev_epj[nj_tot].pos.y = epj[iw][j].pos.y;
	dev_epj[nj_tot].pos.z = epj[iw][j].pos.z;

	dev_epj[nj_tot].vel.x = epj[iw][j].vel.x;
	dev_epj[nj_tot].vel.y = epj[iw][j].vel.y;
	dev_epj[nj_tot].vel.z = epj[iw][j].vel.z;
      
	dev_epj[nj_tot].id() = epj[iw][j].id;
	dev_epj[nj_tot].prop() = epj[iw][j].prop;
      
	nj_tot++;
      }
    }
  }

  template<class RESULT>
  static void CopyToOrigin(const PS::S32 n_walk,
			   const PS::S32 ni[],
			   RESULT* force[]) {
    int n_cnt = 0;
    for(int iw = 0; iw < n_walk; iw++) {
      for(int i = 0; i < ni[iw]; i++) {
	force[iw][i].acc.x = dev_force[n_cnt].acc.x;
	force[iw][i].acc.y = dev_force[n_cnt].acc.y;
	force[iw][i].acc.z = dev_force[n_cnt].acc.z;
	
	force[iw][i].press.x = dev_force[n_cnt].press.x;
	force[iw][i].press.y = dev_force[n_cnt].press.y;
	force[iw][i].press.z = dev_force[n_cnt].press.z;	

	n_cnt++;
      }
    }
  }
};
