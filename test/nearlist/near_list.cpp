#include <iostream>
#include "particle_simulator.hpp"
#include "io_util.hpp"
#include "parameter.hpp"
#include "f_calculator.hpp"
#include <utility>
#include <algorithm>
#include <ctime>

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#error "MPI is not supported yet!"
#endif

constexpr PS::F64 rc = 1.0;
constexpr PS::F64 rc2 = rc * rc;
constexpr PS::S32 num_tot = 3072;
const PS::F64vec box = { 8.0, 8.0, 8.0 };

template<class Tpsys>
void initialize(Tpsys& sys)
{
  PS::MT::init_genrand( static_cast<unsigned long>(time(NULL)));
  //PS::MT::init_genrand(0);
  for(PS::S32 i = 0; i < num_tot; i++) {
    sys[i].id = i;
    sys[i].pos = PS::F64vec(PS::MT::genrand_real1() * box.x,
			    PS::MT::genrand_real1() * box.y,
			    PS::MT::genrand_real1() * box.z);
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
	if(dr2 < rc2) {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	  
	  const PS::S32 tid = omp_get_thread_num();
#else
	  const PS::S32 tid = 0;
#endif
	  nlist[tid].push_back(std::pair<size_t, size_t>(idi, idj) );
	}
      }
    }
  }
};

template<class EPI, class EPJ, class Force>
std::vector<std::vector<std::pair<size_t, size_t> > > CalcForceEpEp_for_nlist<EPI, EPJ, Force>::nlist;

void min_image(PS::F64vec& drij,
	       const PS::F64vec& box)
{
  drij.x -= box.x * std::round(drij.x / box.x);
  drij.y -= box.y * std::round(drij.y / box.y);
  drij.z -= box.z * std::round(drij.z / box.z);
}

template<class Tpsys>
void make_nlist_naive(Tpsys& sys,
		      const PS::F64vec& box,
		      std::vector<std::pair<size_t, size_t> >& nlist_)
{
  PS::S32 num = sys.getNumberOfParticleLocal();
  for(PS::S32 i = 0; i < num; i++) {
    const PS::F64vec ri = sys[i].pos;
    const PS::U32 idi = sys[i].id;
    for(PS::S32 j = 0; j < num; j++) {
      const PS::F64vec rj = sys[j].pos;
      const PS::U32 idj = sys[j].id;
      PS::F64vec drij = ri - rj;
      min_image(drij, box);
      const PS::F64 dr2 = drij * drij;
      if(dr2 < rc2) {
	nlist_.push_back(std::pair<size_t, size_t>(idi, idj) );
      }
    }
  }
}

void compare_nlist(std::vector<std::vector<std::pair<size_t, size_t> > >& nlist,
		   std::vector<std::pair<size_t, size_t> >& nlist_naive)
{
  std::vector<std::pair<size_t, size_t> > merged_list;
  for(size_t i = 0; i < nlist.size(); i++)
    merged_list.insert(merged_list.end(), nlist[i].begin(), nlist[i].end());
  std::sort(merged_list.begin(), merged_list.end());
  std::sort(nlist_naive.begin(), nlist_naive.end());

  if(merged_list.size() != nlist_naive.size()) {
    std::cerr << "nlist_size nlist_naive_size\n";
    std::cerr << merged_list.size() << " " << nlist_naive.size() << std::endl;
    std::cerr << nlist[0].size() << std::endl;
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
  
  std::cout << "Success: Test\n";
}

int main(int argc, char* argv[]) {
  PS::U32 test_num = 3;
  if(argc == 2)
    test_num = atoi(argv[1]);
  else
    std::cerr << "# of test is set to " << test_num << std::endl;

  PS::ParticleSystem<FPDPD> system;
  system.initialize();
  system.setNumberOfParticleLocal(num_tot);  

  PS::DomainInfo dinfo;
  dinfo.initialize(0.3);
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), box);
  dinfo.collectSampleParticle(system);

  PS::TreeForForceShort<RESULT::ForceDPD, EPI::DPD, EPJ::DPD>::Gather tree_prtcl;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
  const int num_th = omp_get_max_threads();
#else
  const int num_th = 1;
#endif
  CalcForceEpEp_for_nlist<EPI::DPD, EPJ::DPD, RESULT::ForceDPD>::nlist.resize(num_th);
  tree_prtcl.initialize(3 * num_tot);
  
  for(PS::U32 i = 0; i < test_num; i++) {
    initialize(system);

    dinfo.decomposeDomain();
    system.exchangeParticle(dinfo);

    tree_prtcl.calcForceAllAndWriteBack(CalcForceEpEp_for_nlist<EPI::DPD, EPJ::DPD, RESULT::ForceDPD>(), system, dinfo);

    std::vector<std::pair<size_t, size_t> > nlist_naive;
    make_nlist_naive(system, box, nlist_naive);
    
    compare_nlist(CalcForceEpEp_for_nlist<EPI::DPD, EPJ::DPD, RESULT::ForceDPD>::nlist,
		  nlist_naive);

    for(int i = 0; i < num_th; i++)
      CalcForceEpEp_for_nlist<EPI::DPD, EPJ::DPD, RESULT::ForceDPD>::nlist[i].clear();
  }  
}
