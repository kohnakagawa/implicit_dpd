#pragma once

#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <stack>
#include <ctime>
#include <cstdlib>
#include <numeric>

// NOTE: We do not support the analysis of constant pressure simulation in this class library.
template<class PtclOrg, class Ptcl>
class PtclConnector {
  // simulation box information
  std::array<double, 3> grid_leng_, box_leng_;
  std::array<int, 3> grid_numb_;
  std::vector<std::vector<int> > ilist_;

  // ptcl data
  std::vector<Ptcl> ptcls_;
  std::vector<PtclOrg> ptcls_org_;
  std::vector<int> ptcl_id_, elem_in_grid_, adrs_grid_;
  
  // for DFS
  std::vector<bool> reg_flag_;
  std::stack<int> id_stack_;
  std::vector<int> patch_id_;

  // near ptcl data
  std::vector<std::vector<int> > near_ptcl_id_;

  // parameters
  double cutof_leng_ = 0.0;

  int GenHash(int* q) {
    if (q[0] < 0 || q[0] > grid_numb_[0] - 1) q[0] -= static_cast<int>(q[0] * grid_numb_[0] / fabs(q[0]));
    if (q[1] < 0 || q[1] > grid_numb_[1] - 1) q[1] -= static_cast<int>(q[1] * grid_numb_[1] / fabs(q[1]));
    if (q[2] < 0 || q[2] > grid_numb_[2] - 1) q[2] -= static_cast<int>(q[2] * grid_numb_[2] / fabs(q[2]));
    return q[0] + grid_numb_[0] * (q[1] + q[2] * grid_numb_[1]);
  }

  template<class array3d>
  void MinImage(array3d& dr) {
    dr[0] -= box_leng_[0] * std::round(dr[0] / box_leng_[0]);
    dr[1] -= box_leng_[1] * std::round(dr[1] / box_leng_[1]);
    dr[2] -= box_leng_[2] * std::round(dr[2] / box_leng_[2]);  
  }
  
  void CreateNearlist() {
    const int all_grid = grid_numb_[0] * grid_numb_[1] * grid_numb_[2];
    ilist_.resize(all_grid);
    for (int i = 0; i < all_grid; i++) {
      ilist_[i].resize(27, 0);
    }
    
    int i_num = 0;
    for (int iz = 0; iz < grid_numb_[2]; iz++) {
      for (int iy = 0; iy < grid_numb_[1]; iy++) {
	for (int ix = 0; ix < grid_numb_[0]; ix++) {
	  int j_num = 0;
	  for (int jz = -1; jz < 2; jz++) {
	    for (int jy = -1; jy < 2; jy++) {
	      for (int jx = -1; jx < 2; jx++) {
		int q[3] = {ix + jx, iy + jy, iz + jz};
		ilist_[i_num][j_num++] = GenHash(q);
	      }
	    }
	  }
	  i_num++;
	}
      }
    }
  }

  int SearchNextPid(const int size) {
    bool flag = true;
    int ret = 0;
    while (flag) {
      int pid = rand() % size;
      if (!reg_flag_[pid]) {
	ret = pid;
	flag = false;
      }
    }
    return ret;
  }

  void DfsWithMod(int pid, const int cur_patch) {
    reg_flag_[pid] = true;
    patch_id_[pid] = cur_patch;
    id_stack_.push(pid);
    while (!id_stack_.empty()) {
      pid = id_stack_.top();
      id_stack_.pop();
      const auto ri = ptcls_[pid].r;
      const auto near_ptcl_size = near_ptcl_id_[pid].size();
      for (auto i = 0u; i < near_ptcl_size; i++) {
	const auto pjd = near_ptcl_id_[pid][i];
	const auto rj = ptcls_[pjd].r;
	auto drij = rj - ri;
	MinImage(drij);
	if (!reg_flag_[pjd]) {
	  ptcls_[pjd].r = ri + drij;
	  reg_flag_[pjd] = true;
	  patch_id_[pjd] = cur_patch;
	  id_stack_.push(pjd);
	}
      }
    }    
  }

  void DfsWithoutMod(int pid, const int cur_patch) {
    reg_flag_[pid] = true;
    patch_id_[pid] = cur_patch;
    id_stack_.push(pid);
    while (!id_stack_.empty()) {
      pid = id_stack_.top();
      id_stack_.pop();
      const auto near_ptcl_size = near_ptcl_id_[pid].size();
      for (auto i = 0u; i < near_ptcl_size; i++) {
	const auto pjd = near_ptcl_id_[pid][i];
	if (!reg_flag_[pjd]) {
	  reg_flag_[pjd] = true;
	  patch_id_[pjd] = cur_patch;
	  id_stack_.push(pjd);
	}
      }
    }
  }
  
  void Dfs(int pid, const int cur_patch, const int mode = 0) {
    switch (mode) {
    case WITH_MOD:
      DfsWithMod(pid, cur_patch);
      break;
    case WITHOUT_MOD:
      DfsWithoutMod(pid, cur_patch);
      break;
    default:
      std::cerr << "Unknown mode.\n";
      std::exit(1);
      break;
    }
  }

public:
  enum {
    WITH_MOD = 0,
    WITHOUT_MOD = 1,
  };

  PtclConnector(const double Lx, const double Ly, const double Lz) {
    box_leng_[0] = Lx;
    box_leng_[1] = Ly;
    box_leng_[2] = Lz;
  }
  
  explicit PtclConnector(const std::array<double, 3>& box_leng) {
    box_leng_ = box_leng;
  }
  
  explicit PtclConnector(const double* box_leng) {
    for (int i = 0; i < 3; i++) {
      box_leng_[i] = box_leng[i];
    }
  }

  void Initialize(const double est_grid_leng,
		  const double cutof_leng) {
    srand((std::size_t)time(nullptr));
    
    assert(est_grid_leng > 0.0);
    assert(cutof_leng > 0.0);
    
    for (int i = 0; i < 3; i++) {
      grid_numb_[i] = static_cast<int>(box_leng_[i] / est_grid_leng);
      grid_leng_[i] = box_leng_[i] / grid_numb_[i];
    }
    
    const int all_grid = grid_numb_[0] * grid_numb_[1] * grid_numb_[2];
    elem_in_grid_.resize(all_grid, 0);
    adrs_grid_.resize(all_grid + 1, 0);
    
    cutof_leng_ = cutof_leng;

    assert(est_grid_leng >= cutof_leng);
    
    CreateNearlist();
  }
  
  bool ReadOneFrame(std::ifstream& fin,
		    std::size_t& time) {
    std::size_t num_ptcl = 0;
    std::string tag, line;
    char buf[4];
    std::getline(fin, line);
    if (fin.eof()) return true;
    num_ptcl = std::stoi(line);
    std::getline(fin, line);
    sscanf(line.c_str(), "%s %u\n", buf, &time);

    ptcls_org_.resize(num_ptcl);
    for (std::size_t i = 0; i < num_ptcl; i++) {
      std::getline(fin, line);
      ptcls_org_[i].readFromString(line.c_str());
    }

    // resize if needed.
    ptcl_id_.resize(num_ptcl, -1);
    patch_id_.resize(num_ptcl, -1);
    near_ptcl_id_.resize(num_ptcl);
  
    // copy to buffer
    ptcls_.clear();
    Ptcl p;
    for (std::size_t i = 0; i < ptcls_org_.size(); i++) {
      p.CopyFromPtclOrg(ptcls_org_[i]);
      if (p.IsTarget()) ptcls_.push_back(p);
    }
    return false;
  }
  
  void RegistPtclIdx() {
    elem_in_grid_.assign(elem_in_grid_.size(), 0);
    adrs_grid_[0] = 0;
    for (std::size_t i = 0; i < ptcls_.size(); i++) {
      ptcls_[i].SetHash(grid_numb_, grid_leng_);
      elem_in_grid_[ptcls_[i].GetHash()]++;
    }

    const int all_grid = grid_numb_[0] * grid_numb_[1] * grid_numb_[2];
  
    for (int i = 0; i < all_grid; i++) {
      adrs_grid_[i + 1] = adrs_grid_[i] + elem_in_grid_[i];
    }
  
    for (std::size_t i = 0; i < ptcls_.size(); i++) {
      const int hash = ptcls_[i].GetHash();
      const int temp_idx = adrs_grid_[hash];
      ptcl_id_[temp_idx] = i;
      adrs_grid_[hash]++;
    }
  
    for (int i = 0; i < all_grid; i++) {
      adrs_grid_[i] = adrs_grid_[i] - elem_in_grid_[i];
    }
  }

  void CheckRegisted() const {
    const int all_grid = grid_numb_[0] * grid_numb_[1] * grid_numb_[2];
    for (int tar = 0; tar < all_grid; tar++) {
      for (int pi = adrs_grid_[tar]; pi < adrs_grid_[tar + 1]; pi++) {
	const int hash = ptcls_[pi].GetHash();
	assert(hash == tar);
      }
    }
  }

  void GenGraph() {
    const std::size_t size = ptcls_.size();
    for (std::size_t i = 0; i < size; i++) {
      const int pi       = ptcl_id_[i];
      const int i_grid   = ptcls_[pi].GetHash();
      const auto ri = ptcls_[pi].r;
      for (int j_grid = 0; j_grid < 27; j_grid++) {
	const int intr_grid = ilist_[i_grid][j_grid];
	const int beg_id    = adrs_grid_[intr_grid    ];
	const int end_id    = adrs_grid_[intr_grid + 1];
	for (int j = beg_id; j < end_id; j++) {
	  const int pj  = ptcl_id_[j];
	  const auto rj = ptcls_[pj].r;
	  auto dr = rj - ri;
	  MinImage(dr);
	  const auto dist = std::sqrt(dr * dr);
	  if ((dist < cutof_leng_) && (pi != pj)) {
	    near_ptcl_id_[pi].push_back(pj);
	  }
	}
      }
    }
  }
  
  // NOTE: O(N^2) method
  void GenGraphNaive() {
    const auto size = ptcls_.size();
    const auto cfl2 = cutof_leng_ * cutof_leng_;
    for (auto i = 0u; i < size; i++) {
      const auto ri = ptcls_[i].r;
      for (auto j = 0u; j < size; j++) {
	const auto rj = ptcls_[j].r;
	auto dr = rj - ri;
	MinImage(dr);
	const auto dist = dr * dr;
	if ((dist < cfl2) and (i != j)) {
	  near_ptcl_id_[i].push_back(j);
	}
      }
    }
  }
  
  int ConnectPtcls(const int mode) {
    const int all_n = ptcls_.size();
    int regist_n = 0, pid = 0, cur_patch = 0;
    const int size  = reg_flag_.size();
    while (regist_n != all_n) {
      pid = SearchNextPid(size);
      Dfs(pid, cur_patch, mode);
      regist_n = std::accumulate(reg_flag_.cbegin(), reg_flag_.cend(), 0);
      std::cout << "regist_n " << regist_n << " ";
      cur_patch++;
    }
    std::cout << std::endl;
    
    return cur_patch;
  }
  
  void ClearForNextStep() {
    ptcl_id_.clear();
    reg_flag_.clear();
    
    for (std::size_t i = 0; i < near_ptcl_id_.size(); i++) {
      near_ptcl_id_[i].clear();
    }
  }

  const std::vector<int>& patch_id() const {
    return patch_id_;
  }

  const std::vector<Ptcl>& ptcls() const {
    return ptcls_;
  }

  std::vector<Ptcl>& ptcls() {
    return ptcls_;
  }

  const std::vector<std::vector<int> >& near_ptcl_id() const {
    return near_ptcl_id_;
  }
};

template<class PtclOrg, class Ptcl>
int gen_connected_image(PtclConnector<PtclOrg, Ptcl>* ptr_connector,
			bool ptcl_data_const = true) {
  ptr_connector->RegistPtclIdx();
  ptr_connector->CheckRegisted();
  
  // generate networks
#ifdef DEBUG
  ptr_connector->GenGraphNaive();
#else
  ptr_connector->GenGraph();
#endif
  
  // carry out DFS
  int num_patch = -1;
  if (ptcl_data_const) {
    num_patch = ptr_connector->ConnectPtcls(PtclConnector<PtclOrg, Ptcl>::WITHOUT_MOD);
  } else {
    num_patch = ptr_connector->ConnectPtcls(PtclConnector<PtclOrg, Ptcl>::WITH_MOD);
  }
  
  return num_patch;
}
