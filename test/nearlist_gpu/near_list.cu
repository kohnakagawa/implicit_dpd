#include <iostream>
#include "particle_simulator.hpp"
#include "io_util.hpp"
#include "parameter.hpp"
#include <utility>
#include <algorithm>
#include <ctime>
#include <iomanip>

constexpr PS::F64 rc = 1.0;
constexpr PS::F64 rc2 = rc * rc;
// const PS::F64vec box = { 32.0, 32.0, 32.0 };
const PS::F64vec box = { 16.0, 16.0, 16.0 };
// const PS::F64vec box = { 8.0, 8.0, 8.0 };


#include "dispatch.hpp"

void min_image(PS::F64vec& drij,
	       const PS::F64vec& box)
{
  drij.x -= box.x * std::round(drij.x / box.x);
  drij.y -= box.y * std::round(drij.y / box.y);
  drij.z -= box.z * std::round(drij.z / box.z);
}

template<class Tpsys>
void initialize(Tpsys& sys)
{
  for(PS::S32 i = 0; i < NUM_TOT; i++) {
    sys[i].id = i;
    sys[i].pos = PS::F64vec(PS::MT::genrand_real1() * box.x,
			    PS::MT::genrand_real1() * box.y,
			    PS::MT::genrand_real1() * box.z);
  }
}

template<typename T>
struct pair_length {
  std::pair<size_t, size_t> pair;
  T length;
};

template<typename T>
struct pred {
  bool operator()(const pair_length<T>& riLeft,
       		  const pair_length<T>& riRight) const {
    return riLeft.pair < riRight.pair;
  }
};

void search_array(const std::vector<pair_length<float> >& gpu_list,
     	          std::vector<pair_length<PS::F64> >& naiv_list)
{
  std::vector<pair_length<PS::F64> > vec;
  for(size_t i = 0; i < naiv_list.size(); i++) {
    std::pair<size_t, size_t> p = naiv_list[i].pair;
    bool flag = false;
    for(size_t j = 0; j < gpu_list.size(); j++) {
      if(p == gpu_list[j].pair)
        flag = true;
    }
    if(!flag)
      vec.push_back(naiv_list[i]);
  }
  std::cout << std::setprecision(15);
  for(size_t i = 0; i < vec.size(); i++) 
    std::cout << vec[i].pair.first << " " << vec[i].pair.second << " " << vec[i].length << std::endl;
}

template<class Tpsys>
void make_nlist_naive(Tpsys& sys,
		      const PS::F64vec& box,
		      std::vector<pair_length<PS::F64> >& nlist_)
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
      if(dr2 < rc2 && idi != idj) {
        pair_length<PS::F64> tmp;
	tmp.pair = std::pair<size_t, size_t>(idi, idj);
	tmp.length = dr2;
	nlist_.push_back(tmp);
      }
    }
  }
}

void compare_nlist(std::vector<pair_length<PS::F64> >& nlist_naive)
{
  std::vector<pair_length<float> > gpu_list;
  
  for(int i = 0; i < NUM_TOT; i++) {
    const uint idi = pairs[i * PAIR_NUM];
    for(int j = 1; j < num_pairs[idi]; j++) {
      pair_length<float> tmp;
      tmp.pair = std::pair<size_t, size_t>(idi, pairs[i * PAIR_NUM + j]);
      tmp.length = length[i * PAIR_NUM + j];
      gpu_list.push_back(tmp);
    }
    for(int j = num_pairs[idi]; j < PAIR_NUM; j++) {
      if(pairs[i * PAIR_NUM + j] != 0xffffffff) {
        std::cout << idi << " " << pairs[i * PAIR_NUM + j] << std::endl;
      }
    }
  }
  
  std::sort(gpu_list.begin(), gpu_list.end(), pred<float>());
  std::sort(nlist_naive.begin(), nlist_naive.end(), pred<PS::F64>());

  std::cout << gpu_list.size() << " " << nlist_naive.size() << std::endl;
  if(gpu_list.size() != nlist_naive.size()) {
    std::cerr << "nlist_size nlist_naive_size\n";
    std::cerr << gpu_list.size() << " " << nlist_naive.size() << std::endl;
    std::cerr << "length of nearlist is different.\n";
    std::cerr << __FILE__ << " " << __LINE__ << std::endl;
        
    search_array(gpu_list, nlist_naive);
        
    PS::Abort();
  }
  
  for(size_t i = 0; i < gpu_list.size(); i++) {
    if(gpu_list[i].pair != nlist_naive[i].pair) {
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

  cudaSetDevice(0);  

  PS::ParticleSystem<FPDPD> system;
  system.initialize();
  system.setNumberOfParticleLocal(NUM_TOT);

  PS::DomainInfo dinfo;
  dinfo.initialize(0.3);
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), box);
  dinfo.collectSampleParticle(system);

  PS::TreeForForceShort<RESULT::ForceDPD, EPI::DPD, EPJ::DPD>::Gather tree_prtcl;
  tree_prtcl.initialize(3 * NUM_TOT);

  //PS::MT::init_genrand( static_cast<unsigned long>(time(NULL)));
  PS::MT::init_genrand(0);
  
  for(PS::U32 i = 0; i < test_num; i++) {
    initialize(system);

    dinfo.decomposeDomain();
    system.exchangeParticle(dinfo);

    //naive method
    std::vector<pair_length<PS::F64> > nlist_naive;
    make_nlist_naive(system, box, nlist_naive);
    
    //
    const PS::S32 n_walk_limit = 200;
    const PS::S32 tag_max = 1;
    
    if(i != 0) {
    	 pairs.set_val(0xffffffff);
    	 length.set_val(0);
    	 num_pairs.set_val(0);
    }
    tree_prtcl.calcForceAllAndWriteBackMultiWalk(DispatchTestKernel,
						 RetrieveTestKernel,
						 tag_max,
						 system,
						 dinfo,
						 n_walk_limit);
    cudaDeviceSynchronize();
    pairs.dev2host(0, NUM_TOT * PAIR_NUM);
    num_pairs.dev2host(0, NUM_TOT);
    length.dev2host(0, NUM_TOT * PAIR_NUM);
    
    compare_nlist(nlist_naive);
  }
  
  pairs.deallocate();
  num_pairs.deallocate();
  ij_disp.deallocate();
  length.deallocate();
  
  Policy::Force::dev_epi.deallocate();
  Policy::Force::dev_epj.deallocate();
  Policy::Force::dev_force.deallocate();
}
