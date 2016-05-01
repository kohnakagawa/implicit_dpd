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
  const char mdstress_errors[10][128] = {
    "the stress type is not correct.",
    "the number of cells at one side is negative.",
    "the local spacing is too small.",
    "the box is not set.",
    "the force decomposition is incorrect.",
    "the number of atoms is not set.",
    "the contribution is not correct",
    "the filename is not set",
    "DistributeInteraction has been called with an incorrect number of atoms",
    "Lapack failed"
  };
  
  PS::F64vec stress_sum_;
  
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
	const auto drji = rj - ri;
	const auto dr2 = drji * drji;
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
      const auto idi   = ptcls_[i].id;
      const auto propi = ptcls_[i].prop;
      const auto densi = ptcls_[i].dens;
      
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
	  const auto inv_dr = 1.0 / dr;
	  const auto propj = ptcls_[j].prop;
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

	  double Force[][3] = {{-0.5 * all_cf * drji.x, -0.5 * all_cf * drji.y, -0.5 * all_cf * drji.z}}; // We do not use action-reaction law.

	  stress_sum_.x += Force[0][0] * drji.x;
	  stress_sum_.y += Force[0][1] * drji.y;
	  stress_sum_.z += Force[0][2] * drji.z;
	  
	  stressgrid.DistributeInteraction(2, atoms_pos, Force, nullptr);
	}
      }
    }
  }

  void AddBondInteractionTerm(mds::StressGrid& stressgrid,
			      const PS::U32 beg_id) {
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
      double Force[][3] = {{-cf_bond * dr01.x, -cf_bond * dr01.y, -cf_bond * dr01.z}};

      stress_sum_.x += Force[0][0] * dr01.x;
      stress_sum_.y += Force[0][1] * dr01.y;
      stress_sum_.z += Force[0][2] * dr01.z;
      
      stressgrid.DistributeInteraction(2, atoms_pos, Force, nullptr);
    }
  }

  void AddAngleInteractionTerm(mds::StressGrid& stressgrid,
			       const PS::U32 beg_id,
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

	const double Ftb0[] = {-cf_bd * (dr12.x - cf_crs[0] * dr01.x),
			       -cf_bd * (dr12.y - cf_crs[0] * dr01.y),
			       -cf_bd * (dr12.z - cf_crs[0] * dr01.z)};
	const double Ftb1[] = {-cf_bd * (dr01.x - cf_crs[1] * dr12.x),
			       -cf_bd * (dr01.y - cf_crs[1] * dr12.y),
			       -cf_bd * (dr01.z - cf_crs[1] * dr12.z)};

	stress_sum_.x += dr01.x * Ftb0[0] + dr12.x * Ftb1[0];
	stress_sum_.y += dr01.y * Ftb0[1] + dr12.y * Ftb1[1];
	stress_sum_.z += dr01.z * Ftb0[2] + dr12.z * Ftb1[2];
	
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

  void AddKineticTerm(mds::StressGrid& stressgrid) {
    const int num_iptcl = ptcls_.size();
    for (int i = 0; i < num_iptcl; i++) {
      double x[]  = {ptcls_[i].r.x, ptcls_[i].r.y, ptcls_[i].r.z};
      double va[] = {ptcls_[i].v.x, ptcls_[i].v.y, ptcls_[i].v.z};
      stressgrid.DistributeKinetic(1.0, x, va, nullptr, -1);
      
      stress_sum_.x -= va[0] * va[0];
      stress_sum_.y -= va[1] * va[1];
      stress_sum_.z -= va[2] * va[2];
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
    stress_sum_.x = stress_sum_.y = stress_sum_.z = 0.0;
    SetDensity(near_ptcl_id);
    AddInterMolecularTerm(stressgrid, near_ptcl_id, nbdf_seed + time);
    AddIntraMolecularTerm(stressgrid, num_mol);
    AddKineticTerm(stressgrid);
    stressgrid.Update();
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
  			  const char* fname,
  			  const double spacing,
  			  const PS::F64vec& box_leng) {
    stressgrid.SetSpacing(spacing);
    stressgrid.SetForceDecomposition(mds_ccfd);
    stressgrid.SetStressType(mds_spat);
    double box[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    for (int i = 0; i < 3; i++) box[i][i] = box_leng[i];
    stressgrid.SetBox(box);
    stressgrid.SetFileName(fname);
    stressgrid.Init();
  }

  void DoAnalysis() override {
    SetSearchRadius(1.2, 1.2);
    ptr_connector->Initialize(est_grid_leng_, cutof_leng_, Parameter::box_leng);

    // for dissipative force calculation
    PS::U32 m_seed = 0;
    const auto fname_seed = cur_dir_ + "/nbdf_rand_seed.txt";
    std::ifstream fin(fname_seed.c_str());
    CHECK_FILE_OPEN(fin, fname_seed);
    fin >> m_seed;
    std::cout << "nbdf_rand_seed is " << m_seed << std::endl;

    // local stress information
    mds::StressGrid stressgrid;
    PS::F64vec ls_box_org(0.0, 0.0, 0.0), ls_box_top(0.0, 0.0, 0.0);

    // sum of press
    const auto fname_psum = cur_dir_ + "/stress_sum.txt";
    std::ofstream fout(fname_psum.c_str());
    
    // main loop
    PS::U32 time = 0;
    AxisAdjuster<Particle> aadjuster;
    bool is_first = true;
    while (true) {
      if (ReadOneFrame(time, [] (const Particle& p) -> bool { return (p.prop == Parameter::Hyphil) || (p.prop == Parameter::Hyphob); })) break;
      std::cout << "time = " << time << std::endl;
      ptr_connector->ResizeIfNeeded(ptcls_.size());
      const auto num_patch = gen_connected_image(ptr_connector.get(), ptcls_, Parameter::box_leng, false);
      std::cout << "# of patch is " << num_patch << std::endl;

#if 0
      // without changing base axis
      if (is_first) {
	const char fname[] = {"loc_stress"}; const double spacing = 5.0;
	SetStressGridParam(stressgrid, fname, spacing, Parameter::box_leng);
	is_first = false;
      }
#else
      // make transform matrix
      const auto& ptch_id = ptr_connector->patch_id();
      aadjuster.CreateMomInertia(ptcls_, ptch_id, GetLargestPatchId(ptch_id, num_patch));
      aadjuster.MakeNewBaseVector(ptcls_);
      aadjuster.ChangeBaseVector(ptcls_);
	
      // calc new box leng
      GetNewBox(ls_box_org, ls_box_top);
	
      // change origin
      ChangeCoordOrigin(ls_box_org);

      if (is_first) {
	const std::string fname =  cur_dir_ + "/loc_stress";
	const double spacing = 5.0;
	std::cout << "ls_box_top: " << ls_box_top << std::endl;
	std::cout << "ls_box_org: " << ls_box_org << std::endl;
	SetStressGridParam(stressgrid, fname.c_str(), spacing, ls_box_top - ls_box_org);
      	is_first = false;
      }
#endif

      CalcLocalStress(stressgrid, ptr_connector->near_ptcl_id(), GetNumOfAmp(), m_seed, time);

      // check mdstress error
      int err_id = 0;
      if ((err_id = stressgrid.GetError()) != 0) {
	std::cerr << "Error occurs in md::stressgrid.\n";
	std::cerr << "Error: " << mdstress_errors[err_id - 1] << std::endl;
	std::exit(1);
      }
      
      fout << std::setprecision(15) << stress_sum_ * Parameter::ibox_leng.x * Parameter::ibox_leng.y * Parameter::ibox_leng.z << std::endl;
      std::cerr << "trace of stress_sum_ is " << stress_sum_.x + stress_sum_.y + stress_sum_.z << std::endl;
      
      ptr_connector->ClearForNextStep();
    }

    DebugDump();

    stressgrid.Write();
  }
};
