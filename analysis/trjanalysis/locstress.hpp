#pragma once
#include "trjanalysis.hpp"
#include "adjust_axis.hpp"
#include "mdstresslib/include/mds_stressgrid.h"
#include "parameter.hpp"
#include "io_util.hpp"
#include "saruprng.hpp"

// static members in Parameter
PS::F64 Parameter::cf_c[Parameter::prop_num][Parameter::prop_num];
PS::F64 Parameter::cf_g[Parameter::prop_num][Parameter::prop_num];
PS::F64 Parameter::cf_r[Parameter::prop_num][Parameter::prop_num];
PS::F64 Parameter::cf_m[Parameter::prop_num][Parameter::prop_num][Parameter::prop_num];
PS::F64 Parameter::cf_s, Parameter::cf_b;

class TrjAnalysisLocStress : public TrjAnalysis<FPDPD, Particle> {
  void SetDensity(const std::vector<std::vector<int> >& near_ptcl_id) {
    const int num_iptcl = ptcls_.size();
    for (int i = 0; i < num_iptcl; i++) {
      const int num_jptcl = near_ptcl_id[i].size();
      const auto ri = ptcls_[i].r;
      double d_sum[Parameter::prop_num] = {0.0};
      
      for (int jj = 0; jj < num_jptcl; jj++) {
	const auto j = near_ptcl_id[i][jj];
	const auto rj = ptcls_[j].r;
	const auto propj = ptcls_[j].prop;
	const auto drij = rj - ri;
	const auto dr2 = drij * drij;
	if (dr2 < Parameter::rc2) {
	  const auto dr = std::sqrt(dr2);
	  d_sum[propj] += (Parameter::rc - dr) * (Parameter::rc - dr);
	}
      }
      
      for (int k = 0; k < Parameter::prop_num; k++)
	ptcls_[i].dens[k] = d_sum[k];
    }
  }

  void AddInterMolecularTerm(mds::StressGrid& stressgrid,
			     const std::vector<std::vector<int> >& near_ptcl_id,
			     const PS::U32 rnd_seed) {
    double atoms_pos[2][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    const int num_iptcl = ptcls_.size();
    for (int i = 0; i < num_iptcl; i++) {
      const int num_jptcl = near_ptcl_id[i].size();
      const auto ri    = ptcls_[i].r;
      const auto vi    = ptcls_[i].v;
      const auto propi = ptcls_[i].prop;
      const auto densi = ptcls_[i].dens;
      const auto idi   = ptcls_[i].id;
      
      for (int j = 0; j < 3; j++) atoms_pos[0][j] = ri[j];
      
      for (int jj = 0; jj < num_jptcl; jj++) {
	const auto j = near_ptcl_id[i][jj];
	const auto rj = ptcls_[j].r;

	const auto idj = ptcls_[j].id;
	const auto drji = rj - ri;
	const auto dr2 = drji * drji;
	if (dr2 < Parameter::rc2) {
	  const auto dvji = ptcls_[j].v - vi;

	  for (int k = 0; k < 3; k++) atoms_pos[1][k] = rj[k];

	  double densij[Parameter::prop_num];
	  for (int k = 0; k < Parameter::prop_num; k++)
	    densij[k] = densi[k] + ptcls_[j].dens[k];

	  auto m_i = idi, m_j = idj;
	  if (idi > idj) {
	    m_i = idj;
	    m_j = idi;
	  }

	  Saru saru(m_i, m_j, rnd_seed);
	  const auto rnd = saru.nrml();

	  const auto dr = std::sqrt(dr2);
	  const auto propj = ptcls_[j].prop;
	  const auto inv_dr = 1.0 / dr;
	  const auto one_m_dr = 1.0 - dr * Parameter::irc;
	  
	  const auto cf_co  = Parameter::cf_c[propi][propj] * (dr - Parameter::arc) * (dr - Parameter::rc) * (dr >= Parameter::arc);
	  double cf_mbd = 0.0;
	  for (PS::S32 k = 0; k < Parameter::prop_num; k++)
	    cf_mbd += densij[k] * Parameter::cf_m[propi][propj][k];
	  cf_mbd *= one_m_dr;

	  const auto wrji = one_m_dr; // pow = 1
	  // const auto wrji = std::sqrt(one_m_dr); // pow = 1 / 2
	  const auto sq_wrji = std::sqrt(wrji);
	  const auto drji_dvji = drji * dvji;
	  const auto all_cf = (cf_co + cf_mbd +  //conservative
			       Parameter::cf_r[propi][propj] * sq_wrji * rnd - //random
			       Parameter::cf_g[propi][propj] * wrji * drji_dvji * inv_dr) * inv_dr; //dissipation

	  double Force[][3] = {{all_cf * drji.x, all_cf * drji.y, all_cf * drji.z}};
	  
	  stressgrid.DistributeInteraction(2, atoms_pos, Force, nullptr);
	}
      }
    }
  }

  void AddBondInteractionTerm(mds::StressGrid& stressgrid,
			      const int beg_id) {
    const auto amp_id = ptcls_[beg_id].amp_id;
    for (auto i = 0u; i < Parameter::all_unit - 1; i++) {
      const auto r0 = ptcls_[beg_id + i    ].r;
      const auto r1 = ptcls_[beg_id + i + 1].r;

      CHECK_EQ(amp_id, ptcls_[beg_id + i + 1].amp_id); // check topology

      const auto dr01 = r1 - r0;
      const auto dr01_2 = dr01 * dr01;
      const auto inv_dr = 1.0 / std::sqrt(dr01_2);
      
      const auto cf_bond = Parameter::cf_spring<Parameter::bond_leng != 0.0>(inv_dr);

      double atoms_pos[][3] = {{r0.x, r0.y, r0.z}, {r1.x, r1.y, r1.z}};
      double Force[][3] = {{cf_bond * dr01.x, cf_bond * dr01.y, cf_bond * dr01.z}};
      
      stressgrid.DistributeInteraction(2, atoms_pos, Force, nullptr);
    }
  }

  void AddAngleInteractionTerm(mds::StressGrid& stressgrid,
			       const int beg_id,
			       const double cf_b) {
    if (Parameter::all_unit > 2) {
      const auto amp_id = ptcls_[beg_id].amp_id;
      for (auto i = 0u; i < Parameter::all_unit - 2; i++) {
	const auto r0 = ptcls_[beg_id + i    ].r;
	const auto r1 = ptcls_[beg_id + i + 1].r;
	const auto r2 = ptcls_[beg_id + i + 2].r;
	
	CHECK_EQ(amp_id, ptcls_[beg_id + i + 1].amp_id);
	CHECK_EQ(amp_id, ptcls_[beg_id + i + 2].amp_id);

	const auto dr01		= r1 - r0;
	const auto dr01_2	= dr01 * dr01;
	const auto inv_dr01	= 1.0 / std::sqrt(dr01_2);

	const auto dr12		= r2 - r1;
	const auto dr12_2	= dr12 * dr12;
	const auto inv_dr12	= 1.0 / std::sqrt(dr12_2);

	const auto inv_dr_prod	= inv_dr01 * inv_dr12;
	const auto cf_bd	= cf_b * inv_dr_prod;
	
	const double inv_dist[]	= {inv_dr01 * inv_dr01, inv_dr12 * inv_dr12};

	const auto in_prod	= dr01 * dr12;
	const double cf_crs[]   = {in_prod * inv_dist[0], in_prod * inv_dist[1]};

	const double Ftb0[] = {cf_bd * (dr12.x - cf_crs[0] * dr01.x),
			       cf_bd * (dr12.y - cf_crs[0] * dr01.y),
			       cf_bd * (dr12.z - cf_crs[0] * dr01.z)};
	const double Ftb1[] = {cf_bd * (dr01.x - cf_crs[1] * dr12.x),
			       cf_bd * (dr01.y - cf_crs[1] * dr12.y),
			       cf_bd * (dr01.z - cf_crs[1] * dr12.z)};
	
	double atoms_pos[][3] = {{r0.x, r0.y, r0.z}, {r1.x, r1.y, r1.z}, {r2.x, r2.y, r2.z}};
	double Force[][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
	
	for (int k = 0; k < 3; k++) {
	  Force[0][k] -= Ftb0[k];
	  Force[1][k] += Ftb0[k] - Ftb1[k];
	  Force[2][k] += Ftb1[k];
	}

	stressgrid.DistributeInteraction(3, atoms_pos, Force, nullptr);
      }
    }
  }


  void SetBendInteractionCoefs(std::vector<double>& cf_bends,
			       const std::vector<PS::U32>& core_amp_id,
			       const PS::U32 num_mol) {
    cf_bends.resize(num_mol);
    const auto num_core_amp_id = core_amp_id.size();
    for (auto i = 0u; i < num_mol; i++) {
      cf_bends[i] = Parameter::cf_b;
      for (auto j = 0u; j < num_core_amp_id; j++) {
	if (i == core_amp_id[j]) {
	  cf_bends[i] = Parameter::cf_b_rigid;
	}
      }
    }
  }

  // ASSUME: ptcls_ data is sorted.
  void AddIntraMolecularTerm(mds::StressGrid& stressgrid,
			     const int num_mol) {
    static std::vector<double> cf_bends;
    static bool is_first = true;
    if (is_first) {
      SetBendInteractionCoefs(cf_bends, ptr_parameter->core_amp_id(), num_mol);
      is_first = false;
    } else {
      // append cf_bends
      const int num_incr = num_mol - cf_bends.size();
      for (int i = 0; i < num_incr; i++)
	cf_bends.push_back(Parameter::cf_b);
    }
    
    for (int i = 0; i < num_mol; i++) {
      AddBondInteractionTerm(stressgrid, Parameter::all_unit * i);
      AddAngleInteractionTerm(stressgrid, Parameter::all_unit * i, cf_bends[i]);
    }
  }

  void CalcLocalStress(mds::StressGrid& stressgrid,
		       const std::vector<std::vector<int> >& near_ptcl_id,
		       const int num_mol,
		       const PS::U32 nbdf_seed,
		       const PS::U32 time) {
    SetDensity(near_ptcl_id);
    AddInterMolecularTerm(stressgrid, near_ptcl_id, nbdf_seed + time);
    AddIntraMolecularTerm(stressgrid, num_mol);
  }

  int GetNumOfAmp() {
    int cnt = 0;
    for (const auto& p : ptcls_) if (p.prop != Parameter::Solvent) cnt++;
    cnt /= Parameter::all_unit;
    return cnt;
  }

  void ShowStressCalcInformation(mds::StressGrid& stressgrid) {
    // For Debug

    //
  }
  
public:
  TrjAnalysisLocStress(const std::string cur_dir,
		       const char* trj_fname,
		       const char* param_fname) : TrjAnalysis(cur_dir, trj_fname, param_fname) {}
  
  ~TrjAnalysisLocStress() override {}

  void SetStressGridParam(mds::StressGrid& stressgrid,
  			  const int grid_n[3],
  			  const char* fname,
  			  const double spacing,
  			  const PS::F64vec& box_leng) {
    stressgrid.SetNumberOfGridCellsX(grid_n[0]);
    stressgrid.SetNumberOfGridCellsY(grid_n[1]);
    stressgrid.SetNumberOfGridCellsZ(grid_n[2]);
    
    stressgrid.SetSpacing(spacing);
    stressgrid.SetForceDecomposition(mds_ccfd);
    stressgrid.SetStressType(mds_spat);
    
    double box[3][3];
    for (int i = 0; i < 3; i++) for(int j = 0; j < 3; j++) box[i][j] = 0.0;
    for (int i = 0; i < 3; i++) box[i][i] = box_leng[i];
    stressgrid.SetBox(box);
    
    stressgrid.SetFileName(fname);

    stressgrid.Init();
  }

  void DoAnalysis() override {
    SetSearchRadius(1.2, 1.2);
    ptr_connector->Initialize(est_grid_leng_, cutof_leng_, Parameter::box_leng);

    // local stress information
    const char fname[] = {"loc_stress"};
    const int grid_n[] = {20, 20, 20};
    const double spacing = 2.0;
    mds::StressGrid stressgrid;
    SetStressGridParam(stressgrid, grid_n, fname, spacing, Parameter::box_leng);

    // for dissipative force calculation
    PS::U32 m_seed = 0;
    const auto fname_seed = cur_dir_ + "/nbdf_rand_seed.txt";
    std::ifstream fin(fname_seed.c_str());
    CHECK_FILE_OPEN(fin, fname_seed);
    fin >> m_seed;
    std::cout << "nbdf_rand_seed is " << m_seed << std::endl;
    
    PS::U32 time = 0;
    AxisAdjuster<Particle> aadjuster;
    while (true) {
      if (ReadOneFrame(time, [] (const Particle& p) -> bool { return (p.prop == Parameter::Hyphil) || (p.prop == Parameter::Hyphob); })) break;
      std::cout << "time = " << time << std::endl;
      ptr_connector->ResizeIfNeeded(ptcls_.size());
      const auto num_patch = gen_connected_image(ptr_connector.get(), ptcls_, Parameter::box_leng, false);
      std::cout << "# of patch is " << num_patch << std::endl;
      
      const auto& ptch_id = ptr_connector->patch_id();
      const auto tar_patch_id = GetLargestPatchId(ptch_id, num_patch);
      
      // change base axis
      aadjuster.CreateMomInertia(ptcls_, ptch_id, tar_patch_id);
      aadjuster.DoTransform(ptcls_);
      
      CalcLocalStress(stressgrid, ptr_connector->near_ptcl_id(), GetNumOfAmp(), m_seed, time);
      
      // DebugDump(ptcls_, time);

      ptr_connector->ClearForNextStep();
    }
  }
};