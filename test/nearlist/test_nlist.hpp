#pragma once

#include <utility>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>

namespace test_nlist {
  static constexpr PS::F64 rc = 1.0;
  static constexpr PS::F64 rc2 = rc * rc;
  static constexpr PS::S32 num_tot = 5184;
  static constexpr PS::F64vec box(12.0, 12.0, 12.0);
}

namespace TN = test_nlist;

template<class Tpsys>
void initialize(Tpsys& sys,
		PS::DomainInfo& dinfo)
{
  sys.initialize();
  sys.setNumberOfParticleLocal(TN::num_tot);
  
  dinfo.initialize(0.3);
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), TN::box);
  dinfo.collectSampleParticle(sys);
  dinfo.decomposeDomain();
  sys.exchangeParticle(dinfo);
  
  for(PS::S32 i = 0; i < TN::num_tot; i++) {
    sys[i].id = i;
    sys[i].pos = PS::FP64vec(PS::MT::genrand_real1() * TN::box.x,
			     PS::MT::genrand_real1() * TN::box.y,
			     PS::MT::genrand_real1() * TN::box.z);
  }
}

template<class EPI, class EPJ, class Force>
struct CalcForceEpEp_for_nlist {
  static std::vector<std::vector<std::pair<size_t, size_t> > > nlist;
  void operator() (const EPI *epi,
		   const PS::S32 ni,
		   const EPJ *epj,
		   const PS::S32 nj,
		   Force * result)
  {
    for(PS::S32 i = 0; i < ni; i++) {
      const PS::F64vec ri = epi[i].pos;
      const PS::U32 idi = epi[i].id;
      for(PS::S32 j = 0; j < nj; j++) {
	const PS::F64vec rj = epj[j].pos;
	const PS::U32 idj = epj[j].id;
	const PS::F64vec drij = ri - rj;
	const PS::F64 dr2 = drij * drij;
	if(dr2 < TN::rc2) {
	  const int tid = omp_get_thread_num();
	  nlist[tid].push_back(std::pair<size_t, size_t>(idi, idj) );
	}
      }
    }
  }
};

template<class Tpsys, class EPI, class EPJ, class Force>
void make_nlist_FDPS(Tpsys& sys, PS::DomainInfo& dinfo) {
  const int num_tid = omp_get_max_threads();
  CalcForceEpEp_for_nlist<EPI, EPJ, Force>::nlist.resize(num_tid);
  PS::TreeForForceShort<Force, EPI, EPJ>::Scatter tree;
  tree.initialize(3 * TN::num_tot);
  tree.calcForceAllAndWriteBack(CalcForceEpEp_for_nlist(), sys, dinfo);
}

void min_image(PS::F64vec& drij,
	       const PS::F64vec& box)
{
  drij.x -= box.x * std::round(drij.x / box.x);
  drij.y -= box.y * std::round(drij.y / box.y);
  drij.z -= box.z * std::round(drij.z / box.z);
}

template<class Tpsys>
void make_nlist_naive(Tpsys& sys,
		      PS::F64vec& box,
		      std::vector<std::pair<size_t, size_t> >& nlist) 
{
  PS::S32 num = sys.getNumberOfParticleLocao();
  for(PS::S32 i = 0; i < num; i++) {
    const PS::F64vec ri = sys[i].pos;
    const PS::U32 idi = sys[i].id;
    for(PS::S32 j = i + 1; j < num; j++) {
      const PS::F64vec rj = sys[i].pos;
      const PS::U32 idj = sys[i].id;
      PS::F64vec drij = ri - rj; min_image(drij, box);
      PS::F64 dr2 = drij * drij;
      if(dr2 < TN::rc2) {
	nlist.push_back(std::pair<size_t, size_t>(idi, idj) );
      }
    }
  }
}

void compare_nlist(const std::vector<std::vector<std::pair<size_t, size_t> > >& nlist,
		   const std::vector<std::pair<size_t, size_t> >& nlist_naive)
{
  std::vector<std::pair<size_t, size_t> > merged_list;
  for(size_t i = 0; i < nlist.size(); i++)
    merged_list.insert(merged_list.end(), nlist[i].begin(), nlist[i].end());
  std::sort(merged_list.begin(), merged_list.end());
  std::sort(nlist_naive.begin(), nlist_naive.end());

  if(nlist.size() != nlist_naive.size()) {
    std::cerr << "length of nearlist is different.\n";
    std::cerr << __FILE__ << " " << __LINE__ << std::endl;
    PS::Abort();
  }
  
  for(size_t i = 0; i < nlist.size(); i++) {
    if(merged_list[i] != nlist_naive[i]) {
      std::cerr << "pair list is different\n";
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;
      PS::Abort();
    }
  }
  
  std::cout << "Success: nearlist test\n";
}

void run_check() {
  PS::ParticleSystem<FPDPD> system;
  PS::DomainInfo dinfo;
  initialize(system, dinfo);
  make_nlist_FDPS<FPDPD, EPIDPD, EPJDPD, ForceDPD>(system, dinfo);
  std::vector<std::pair<size_t, size_t> > nlist_naive;
  make_nlist_naive(system, TN::box, nlist_naive);
  compare_nlist(CalcForceEpEp_for_nlist<EPIDPD, EPJDPD, ForceDPD>::nlist,
		nlist_naive);
}
