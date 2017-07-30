#include <iostream>
#include "particle_simulator.hpp"
#include "io_util.hpp"
#include "parameter.hpp"
#include "f_calculator.hpp"
#include <utility>
#include <algorithm>
#include <ctime>

constexpr char Parameter::atom_type[21];

constexpr PS::F64 rc = 1.0;
constexpr PS::F64 rc2 = rc * rc;
constexpr PS::S32 ini_num_tot = 1536;
constexpr PS::S32 increase_rate = 20;
const PS::F64vec box = { 8.0, 8.0, 8.0 };

template<class Tpsys>
void initialize(Tpsys& sys)
{
  PS::MT::init_genrand(10);
  for(PS::S32 i = 0; i < ini_num_tot; i++) {
    sys[i].id = i;
    sys[i].prop = sys[i].amp_id = sys[i].unit = 0;
    sys[i].pos = PS::F64vec(PS::MT::genrand_real1() * box.x,
                            PS::MT::genrand_real1() * box.y,
                            PS::MT::genrand_real1() * box.z);
  }
}

template<class Tpsys>
void increase_particles(Tpsys& sys)
{
  const PS::S32 num_proc = PS::Comm::getNumberOfProc();
  const PS::S32 rank     = PS::Comm::getRank();
  static PS::S32 beg = (increase_rate / num_proc) * rank + ini_num_tot;
  static PS::S32 end = (increase_rate / num_proc) * (rank + 1) + ini_num_tot;

  FPDPD new_ptcl;
  for(PS::S32 i = beg; i < end; i++) {
    new_ptcl.id = i;
    sys[i].prop = sys[i].amp_id = sys[i].unit = 0;
    new_ptcl.pos = PS::F64vec((i - ini_num_tot) * 0.0273,
                              (i - ini_num_tot) * 0.0232,
                              (i - ini_num_tot) * 0.0211);
    if (new_ptcl.pos.x > box.x) new_ptcl.pos.x -= box.x;
    if (new_ptcl.pos.y > box.y) new_ptcl.pos.y -= box.y;
    if (new_ptcl.pos.z > box.z) new_ptcl.pos.z -= box.z;

    sys.AddNewParticle(new_ptcl);
  }

  beg += increase_rate;
  end += increase_rate;
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
      if (dr2 < rc2) {
        nlist_.push_back(std::pair<size_t, size_t>(idi, idj) );
      }
    }
  }
}

void compare_nlist(std::vector<std::vector<std::pair<size_t, size_t> > >& nlist,
                   std::vector<std::pair<size_t, size_t> >& nlist_naive,
                   std::vector<std::pair<size_t, size_t> >& merged_list_loc)
{
  for(size_t i = 0; i < nlist.size(); i++)
    merged_list_loc.insert(merged_list_loc.end(), nlist[i].begin(), nlist[i].end());
  std::sort(merged_list_loc.begin(), merged_list_loc.end());
  std::sort(nlist_naive.begin(), nlist_naive.end());

#ifndef PARTICLE_SIMULATOR_MPI_PARALLEL
  if(merged_list_loc.size() != nlist_naive.size()) {
    std::cerr << "nlist_size nlist_naive_size\n";
    std::cerr << merged_list_loc.size() << " " << nlist_naive.size() << std::endl;
    std::cerr << nlist[0].size() << std::endl;
    std::cerr << "length of nearlist is different.\n";
    std::cerr << __FILE__ << " " << __LINE__ << std::endl;
    PS::Abort();
  }

  for(size_t i = 0; i < nlist.size(); i++) {
    if(merged_list_loc[i] != nlist_naive[i]) {
      std::cerr << "pair list is different\n";
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;
      PS::Abort();
    }
  }

  std::cout << "Success: Test\n";
#endif
}

void show_merged_list_glb(std::vector<std::pair<size_t, size_t>>& merged_list_loc,
                          std::ofstream& fout) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  const PS::S32 num_proc     = PS::Comm::getNumberOfProc();
  const PS::S32 len_list_loc = merged_list_loc.size();
  const PS::S32 len_list_glb = PS::Comm::getSum(len_list_loc);

  std::vector<std::pair<size_t, size_t> > merged_list_glb(len_list_glb);
  std::vector<int> n_list(num_proc);
  std::vector<int> n_list_disp(num_proc + 1);

  PS::Comm::allGather(&len_list_loc, 1, n_list.data());
  n_list_disp[0] = 0;
  for (PS::S32 i = 0; i < num_proc; i++) {
    n_list_disp[i + 1] = n_list_disp[i] + n_list[i];
  }

  PS::Comm::allGatherV(merged_list_loc.data(), len_list_loc,
                       merged_list_glb.data(), n_list.data(), n_list_disp.data());

  if (PS::Comm::getRank() == 0) {
    std::sort(merged_list_glb.begin(), merged_list_glb.end());
    for (PS::S32 i = 0; i < len_list_glb; i++) {
      fout << merged_list_glb[i].first
           << " "
           << merged_list_glb[i].second
           << std::endl;
    }
  }
#else
  const PS::S32 len = merged_list_loc.size();
  std::sort(merged_list_loc.begin(), merged_list_loc.end());
  for (PS::S32 i = 0; i < len; i++) {
    fout << merged_list_loc[i].first
         << " "
         << merged_list_loc[i].second
         << std::endl;
  }
#endif
}

int main(int argc, char* argv[]) {
  PS::Initialize(argc, argv);

  PS::U32 test_num = 20;
  if (argc == 2) {
    test_num = std::atoi(argv[1]);
  } else {
    std::cerr << "# of test is set to " << test_num << std::endl;
  }

  PS::ParticleSystem<FPDPD> system;
  system.initialize();
  if (PS::Comm::getRank() == 0) {
    system.setNumberOfParticleLocal(ini_num_tot);
    initialize(system);
  } else {
    system.setNumberOfParticleLocal(0);
  }

  PS::DomainInfo dinfo;
  dinfo.initialize(0.3);
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), box);
  dinfo.decomposeDomainAll(system);
  system.exchangeParticle(dinfo);

  PS::TreeForForceShort<RESULT::ForceDPD, EPI::DPD, EPJ::DPD>::Gather tree_prtcl;
  const int num_th = PS::Comm::getNumberOfThread();

  CalcForceEpEp_for_nlist<EPI::DPD, EPJ::DPD, RESULT::ForceDPD>::nlist.resize(num_th);
  tree_prtcl.initialize(10 * ini_num_tot);

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  std::ofstream fout("list_mpi.dat");
#else
  std::ofstream fout("list.dat");
#endif

  for (PS::U32 i = 0; i < test_num; i++) {
    std::cout << i << std::endl;

    //increase # of particles
    increase_particles(system);

    dinfo.collectSampleParticle(system);

    dinfo.decomposeDomainAll(system);
    system.exchangeParticle(dinfo);

    tree_prtcl.calcForceAllAndWriteBack(CalcForceEpEp_for_nlist<EPI::DPD, EPJ::DPD, RESULT::ForceDPD>(), system, dinfo);

    std::vector<std::pair<size_t, size_t> > nlist_naive;
#ifndef PARTICLE_SIMULATOR_MPI_PARALLEL
    make_nlist_naive(system, box, nlist_naive);
#endif

    std::vector<std::pair<size_t, size_t> > merged_list_loc; // merged_list for each proc
    compare_nlist(CalcForceEpEp_for_nlist<EPI::DPD, EPJ::DPD, RESULT::ForceDPD>::nlist,
                  nlist_naive,
                  merged_list_loc);

    show_merged_list_glb(merged_list_loc, fout);

    for (int i = 0; i < num_th; i++) {
      CalcForceEpEp_for_nlist<EPI::DPD, EPJ::DPD, RESULT::ForceDPD>::nlist[i].clear();
    }
  }

  PS::Finalize();
}
