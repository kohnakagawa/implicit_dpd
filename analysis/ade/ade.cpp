#include <iostream>
#include <cassert>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <stack>
#include <array>
#include <limits>
#include "particle_simulator.hpp"
#include "io_util.hpp"
#include "parameter.hpp"
#include "ptcl_class.hpp"

#define CHECK_FILE_OPEN(fin, name)			\
  do {							\
    if (!fin) {						\
      cerr << "Cannot open file " << name << std::endl;	\
      cerr << __FILE__ << " " << __LINE__ << std::endl;	\
      exit(1);						\
    }							\
  } while (false)

#define CHECK_EQ(val0, val1)						\
  do {									\
    if (val0 != val1) {							\
      std::cerr << "Is not equal.\n";					\
      std::cerr << "at " __FILE__ << " " << __LINE__ << std::endl;	\
      std::cerr << #val0 << " " << #val1 << std::endl;			\
      std::cerr << val0 << " " << val1 << std::endl;			\
      std::exit(1);							\
    }									\
  } while (false)

#define CHECK_LE(val0, val1)						\
  do {									\
    if (val0 > val1) {							\
      std::cerr << "val0 <= val1 is not statisfied.\n";			\
      std::cerr << "at " __FILE__ << " " << __LINE__ << std::endl;	\
      std::cerr << #val0 << " " << #val1 << std::endl;			\
      std::cerr << val0 << " " << val1 << std::endl;			\
      std::exit(1);							\
    }									\
  } while (false)							\

using namespace std;

constexpr double cutof_leng = 0.9; // this value shoule be less than membrane thickness.

struct Particle {
  PS::F64vec r;
  int prop, amp_id, unit, ptcl_in_cell;
};

void set_scale(PS::F64vec& box_leng,
	       array<int, 3>& grid_numb,
	       PS::F64vec& grid_leng) {
  for (int i = 0; i < 3; i++) {
    grid_numb[i] = static_cast<int>(box_leng[i] / cutof_leng);
    grid_leng[i] = box_leng[i] / grid_numb[i];
  }
}

int gen_hash(const PS::F64vec& r,
	     const array<int, 3>& grid_numb,
	     const PS::F64vec& grid_leng) {
  int i = static_cast<int>(r.x / grid_leng.x), j = static_cast<int>(r.y / grid_leng.y), k = static_cast<int>(r.z / grid_leng.z);
 
  if(i == grid_numb[0]) i--;
  if(j == grid_numb[1]) j--;
  if(k == grid_numb[2]) k--;
  
  const int hash = i + grid_numb[0] * (j + k * grid_numb[1]);
  assert((hash >= 0 && hash < grid_numb[0] * grid_numb[1] * grid_numb[2]) && "hash value is out of range.");
  return hash;
}

int gen_hash(int* q,
	     const array<int, 3>& grid_numb) {
  if(q[0] < 0 || q[0] > grid_numb[0] - 1) q[0] -= static_cast<int>(q[0] * grid_numb[0] / fabs(q[0]));
  if(q[1] < 0 || q[1] > grid_numb[1] - 1) q[1] -= static_cast<int>(q[1] * grid_numb[1] / fabs(q[1]));
  if(q[2] < 0 || q[2] > grid_numb[2] - 1) q[2] -= static_cast<int>(q[2] * grid_numb[2] / fabs(q[2]));
  return q[0] + grid_numb[0] * (q[1] + q[2] * grid_numb[1]);
}

void min_image(PS::F64vec& dr,
	       const PS::F64vec& box) {
  // dr.x -= box.x * round(dr.x / box.x);
  // dr.y -= box.y * round(dr.y / box.y);
  // dr.z -= box.z * round(dr.z / box.z);  
}

void create_nearlist(vector<vector<int> >& ilist,
		     const array<int, 3>& grid_numb) {
  const int all_grid = ilist.size();
  for (int i = 0; i < all_grid; i++) {
    ilist[i].resize(27);
  }
  
  int i_num = 0;
  for (int iz = 0; iz < grid_numb[2]; iz++) {
    for (int iy = 0; iy < grid_numb[1]; iy++) {
      for (int ix = 0; ix < grid_numb[0]; ix++) {
	int j_num = 0;
	for (int jz = -1; jz < 2; jz++) {
	  for (int jy = -1; jy < 2; jy++) {
	    for (int jx = -1; jx < 2; jx++) {
	      int q[3] = {ix + jx, iy + jy, iz + jz};
	      ilist[i_num][j_num++] = gen_hash(q, grid_numb);
	    }
	  }
	}
	i_num++;
      }
    }
  }
}

void regist_ptcl_idx(vector<int>& ptcl_id,
		     vector<int>& elem_in_grid,
		     vector<int>& adrs_grid,
		     vector<Particle>& ptcls,
		     const PS::F64vec& grid_leng,
		     const PS::F64vec& box_leng,
		     const array<int, 3>& grid_numb) {
  elem_in_grid.assign(elem_in_grid.size(), 0);
  adrs_grid[0] = 0;
  for (size_t i = 0; i < ptcls.size(); i++) {
    const int hash = gen_hash(ptcls[i].r, grid_numb, grid_leng);
    ptcls[i].ptcl_in_cell = hash;
    elem_in_grid[hash]++;
  }

  const int all_grid = grid_numb[0] * grid_numb[1] * grid_numb[2];
  
  for (int i = 0; i < all_grid; i++) {
    adrs_grid[i + 1] = adrs_grid[i] + elem_in_grid[i];
  }
  
  for (size_t i = 0; i < ptcls.size(); i++) {
    const int hash = ptcls[i].ptcl_in_cell;
    const int temp_idx = adrs_grid[hash];
    ptcl_id[temp_idx] = i;
    adrs_grid[hash]++;
  }
  
  for (int i = 0; i < all_grid; i++) {
    adrs_grid[i] = adrs_grid[i] - elem_in_grid[i];
  }
}

void check_reg(const vector<int>& adrs_grid, 
	       const vector<int>& ptcl_id,
	       const vector<Particle>& ptcls,
	       const int all_grid,
	       const array<int, 3>& grid_numb,
	       const PS::F64vec& grid_leng)
{
  for (int tar = 0; tar < all_grid; tar++) {
    for (int pi = adrs_grid[tar]; pi < adrs_grid[tar + 1]; pi++) {
      const int hash = gen_hash(ptcls[ptcl_id[pi]].r, grid_numb, grid_leng);
      assert(hash == tar);
    }
  }
}

void generate_graph(vector<Particle>& ptcls,
		    vector<int>& ptcl_id,
		    vector<int>& adrs_grid,
		    vector<vector<int> >& near_ptcl_id,
		    vector<vector<int> >& ilist,
		    const PS::F64vec& box_leng) {
  const size_t size = ptcls.size();
  for (size_t i = 0; i < size; i++) {
    const int pi       = ptcl_id[i];
    const int i_grid   = ptcls[pi].ptcl_in_cell;
    const auto ri = ptcls[pi].r;
    for (int j_grid = 0; j_grid < 27; j_grid++) {
      const int intr_grid = ilist[i_grid][j_grid];
      const int beg_id    = adrs_grid[intr_grid    ];
      const int end_id    = adrs_grid[intr_grid + 1];
      for (int j = beg_id; j < end_id; j++) {
	const int pj  = ptcl_id[j];
	const auto rj = ptcls[pj].r;
	auto dr = rj - ri;
	min_image(dr, box_leng);
	const auto dist = std::sqrt(dr * dr);
	if ((dist < cutof_leng) && (pi != pj)) {
	  near_ptcl_id[pi].push_back(pj);
	}
      }
    }
  }
}

void generate_graph_naive(vector<Particle>& ptcls,
			  vector<vector<int> >& near_ptcl_id,
			  const PS::F64vec& box_leng) {
  const auto size = ptcls.size();
  const auto cfl2 = cutof_leng * cutof_leng;
  for (auto i = 0u; i < size; i++) {
    const auto ri = ptcls[i].r;
    for (auto j = 0u; j < size; j++) {
      const auto rj = ptcls[j].r;
      auto dr = rj - ri;
      min_image(dr, box_leng);
      const auto dist = dr * dr;
      if ((dist < cfl2) and (i != j)) {
	near_ptcl_id[i].push_back(j);
      }
    }
  }
}

void check(vector<vector<int> >& cell,
	   vector<vector<int> >& naive) {
  std::vector<std::pair<int, int> > nlist_cell;
  std::vector<std::pair<int, int> > nlist_naive;
  
  for (auto i = 0u; i < cell.size(); i++)
    for (auto j = 0u; j < cell[i].size(); j++)
      nlist_cell.push_back({i, cell[i][j]});

  for (auto i = 0u; i < naive.size(); i++)
    for (auto j = 0u; j < naive[i].size(); j++)
      nlist_naive.push_back({i, naive[i][j]});

  std::sort(nlist_cell.begin(), nlist_cell.end());
  std::sort(nlist_naive.begin(), nlist_naive.end());

  CHECK_EQ(nlist_cell.size(), nlist_naive.size());
  for (auto i = 0u; i < nlist_naive.size(); i++)
    assert(nlist_cell[i] == nlist_naive[i]);
}

void check_nearlist_correct(const vector<Particle>& ptcls,
			    const vector<vector<int> >& near_ptcl_id,
			    const PS::F64vec& box_leng) {
  const auto size = ptcls.size();
  const auto cfl2 = cutof_leng * cutof_leng;
  for (auto i = 0u; i < size; i++) {
    const auto ri = ptcls[i].r;
    for (auto jj = 0u; jj < near_ptcl_id[i].size(); jj++) {
      const auto j = near_ptcl_id[i][jj];
      const auto rj = ptcls[j].r;
      auto dr = rj - ri;
      min_image(dr, box_leng);
      const auto dist = dr * dr;
      CHECK_LE(dist, cfl2);
    }
  }
}

void dfs(int pid,
	 vector<bool>& reg_flag,
	 const vector<vector<int> >& near_ptcl_id,
	 vector<int>& patch_id,
	 vector<Particle>& ptcls,
	 stack<int>& id_stack,
	 const int cur_patch) {
  reg_flag[pid] = true;
  patch_id[pid] = cur_patch;
  id_stack.push(pid);
  while (!id_stack.empty()) {
    pid = id_stack.top();
    id_stack.pop();
    const auto ri = ptcls[pid].r;
    const auto near_ptcl_size = near_ptcl_id[pid].size();
    for (auto i = 0u; i < near_ptcl_size; i++) {
      const auto pjd = near_ptcl_id[pid][i];
      const auto rj = ptcls[pjd].r;
      const auto drij = ri - rj;
      const auto dr_norm = std::sqrt(drij * drij);
      CHECK_LE(dr_norm, cutof_leng);
      if (!reg_flag[pjd]) {
	reg_flag[pjd] = true;
	patch_id[pjd] = cur_patch;
	id_stack.push(pjd);
      }
    }
  }
}

int search_next_pid(vector<bool>& reg_flag,
		    const int size,
		    vector<Particle>& ptcls) {
  bool flag = true;
  int ret = 0;
  while (flag) {
    int pid = rand() % size;
    if (!reg_flag[pid]) {
      ret = pid;
      flag = false;
    }
  }
  return ret;
}

void search_patch(vector<bool>& reg_flag,
		  const vector<vector<int> >& near_ptcl_id,
		  vector<int>& patch_id,
		  vector<Particle>& ptcls,
		  stack<int>& id_stack,
		  int& cur_patch) {
  const int all_n = ptcls.size();
  int    regist_n = 0;
  int    pid      = 0;
  const int size  = reg_flag.size();
  while (regist_n != all_n) {
    pid = search_next_pid(reg_flag, size, ptcls);
    dfs(pid, reg_flag, near_ptcl_id, patch_id, ptcls, id_stack, cur_patch);
    regist_n = std::accumulate(reg_flag.cbegin(), reg_flag.cend(), 0);
    cout << "regist_n " << regist_n << " ";
    cur_patch++;
  }
  cout << endl;
}

void print_inout_elem(const vector<int>& patch_id,
		      ofstream& fout,
		      const int cur_patch,
		      const int cur_time,
		      const vector<Particle>& ptcls) {
  vector<int> ret;
  for (int i = 0; i < cur_patch; i++) {
    ret.push_back(count(patch_id.begin(), patch_id.end(), i));
  }

  sort(ret.begin(), ret.end(), greater<int>());
  if (ret.size() == 2) {
    fout << cur_time << " " << ret[1] << " " << ret[0] << endl;
  }
}

bool read_one_frame(vector<FPDPD>& ptcls_org,
		    vector<Particle>& ptcls,
		    std::ifstream& fin,
		    size_t& time) {
  size_t num_ptcl = 0;
  string tag, line;
  char buf[4];
  std::getline(fin, line);
  if (fin.eof()) return true;
  num_ptcl = std::stoi(line);
  std::getline(fin, line);
  sscanf(line.c_str(), "%s %u\n", buf, &time);
  
  ptcls_org.resize(num_ptcl);
  for (size_t i = 0; i < num_ptcl; i++) {
    std::getline(fin, line);
    FPDPD pt;
    sscanf(line.c_str(), "%c %lf %lf %lf %u %u %u %u %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	   &(buf[0]), &(pt.pos.x), &(pt.pos.y), &(pt.pos.z),
	   &(pt.id), &(pt.prop), &(pt.amp_id), &(pt.unit),
	   &(pt.vel.x), &(pt.vel.y), &(pt.vel.z), &(pt.vel_buf.x), &(pt.vel_buf.y), &(pt.vel_buf.z),
	   &(pt.acc.x), &(pt.acc.y), &(pt.acc.z));
    ptcls_org[i] = pt;
  }
  
  // copy to buffer
  ptcls.clear();
  for (size_t i = 0; i < ptcls_org.size(); i++) {
    Particle p;
    p.r			= ptcls_org[i].pos;
    p.prop		= ptcls_org[i].prop;
    p.amp_id		= ptcls_org[i].amp_id;
    p.unit		= ptcls_org[i].unit;
    p.ptcl_in_cell	= -1;
    if (p.unit == 0) {
      ptcls.push_back(p);
    }
  }
  return false;
}

void dump_one_frame(const vector<Particle>& ptcls,
		    const vector<int>& patch_id) {
  std::ofstream fout("one_frame.xyz");
  fout << ptcls.size() << std::endl;;
  fout << "one frame\n";
  for (auto i = 0u; i < ptcls.size(); i++) {
    if (patch_id[i] == 0) {
      fout << "O " << ptcls[i].r << " " << i << std::endl;
    } else {
      fout << "N " << ptcls[i].r << " " << i << std::endl;
    }
  }
}

void dump_radial_hist(const vector<Particle>& ptcls) {
  ofstream fout("radial_dist.txt");
  PS::F64vec cm_pos(0.0, 0.0, 0.0);
  for (const auto& p : ptcls) {
    cm_pos += p.r;
  }
  cm_pos /= ptcls.size();
  
  for (const auto& p : ptcls) {
    const auto cm2ptcl = p.r - cm_pos;
    const auto dist = std::sqrt(cm2ptcl * cm2ptcl);
    fout << dist << std::endl;
  }
}

void dump_near_ptcl(const vector<Particle>& ptcls,
		    const vector<vector<int> >& near_ptcl_id) {
  ofstream fout("near_ptcl_config.xyz");
  const int target_id = 33590;
  
  fout << 1 + near_ptcl_id[target_id].size() << std::endl;
  fout << "near ptcls\n";
  fout << "C " << ptcls[target_id].r << std::endl;
  for (size_t i = 0; i < near_ptcl_id[target_id].size(); i++) {
    fout << "O " << ptcls[near_ptcl_id[target_id][i]].r << std::endl;
  }
}

void calc_min_hei(const vector<Particle>& ptcls,
		  const vector<int>& patch_id) {
  vector<Particle> inout[2];
  for (size_t i = 0; i < ptcls.size(); i++) {
    if (patch_id[i] == 0) inout[0].push_back(ptcls[i]);
    if (patch_id[i] == 1) inout[1].push_back(ptcls[i]);
  }
  
  ofstream fout("min_height.txt");
  for (size_t i = 0; i < inout[0].size(); i++) {
    const auto ri = inout[0][i].r;
    double min_hei = numeric_limits<double>::max();
    for (size_t j = 0; j < inout[1].size(); j++) {
      const auto rj = inout[1][j].r;
      const auto dr = rj - ri;
      const double dr_norm = std::sqrt(dr * dr);
      if (dr_norm < min_hei) {
	min_hei = dr_norm;
      }
    }
    fout << min_hei << std::endl;
  }
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
    cerr << "argv[1] == target directory \n";
    exit(1);
  }

  srand((unsigned int)time(nullptr));
  const string s_cur  = argv[1];
  const string fname[] = {
    s_cur + "/traject.xyz",
    s_cur + "/run_param.txt",
    s_cur + "/area_elem.txt",
  };
  
  // open traj data
  vector<FPDPD> ptcls_org;
  vector<Particle> ptcls;
  size_t time = 0;
  ifstream fin_p(fname[0].c_str());
  CHECK_FILE_OPEN(fin_p, fname[0]);

  // read parameters & set scales
  ifstream f_m(fname[1].c_str());
  CHECK_FILE_OPEN(f_m, fname[1]);
  map<string, vector<string> > tag_val;
  Parameter::ReadTagValues(f_m, tag_val);
  PS::F64vec box_leng, grid_leng;
  Parameter::Matching(&(box_leng[0]), std::string("box_leng"), tag_val, 3);
  array<int, 3> grid_numb;
  set_scale(box_leng, grid_numb, grid_leng);
  std::cerr << "Read parameters.\n";

  // make nearlist
  const int all_grid = grid_numb[0] * grid_numb[1] * grid_numb[2];
  vector<int> elem_in_grid, adrs_grid;
  elem_in_grid.resize(all_grid, 0);
  adrs_grid.resize(all_grid + 1, 0);
  vector<vector<int> > ilist(all_grid);
  create_nearlist(ilist, grid_numb);
  std::cerr << "Make nealist.\n";

  // main loop
  ofstream fout(fname[2].c_str());
  vector<int> ptcl_id, patch_id;
  vector<vector<int> > near_ptcl_id, near_ptcl_id_naive;
  vector<bool> reg_flag;
  stack<int> id_stack;
  while (true) {
    if (read_one_frame(ptcls_org, ptcls, fin_p, time)) break;
    std::cout << "time = " << time << std::endl;
    cerr << "Read particles data.\n";
    cerr << ptcls.size() << " " << ptcls_org.size() << endl;
    
    const size_t p_size = ptcls.size();
    ptcl_id.resize(p_size, -1);
    patch_id.resize(p_size, -1);
    near_ptcl_id.resize(p_size);
    near_ptcl_id_naive.resize(p_size);
    
    regist_ptcl_idx(ptcl_id, elem_in_grid, adrs_grid, ptcls, grid_leng, box_leng, grid_numb);
    check_reg(adrs_grid, ptcl_id, ptcls, all_grid, grid_numb, grid_leng);
      
    // generate_graph(ptcls, ptcl_id, adrs_grid, near_ptcl_id, ilist, box_leng);
    generate_graph_naive(ptcls, near_ptcl_id, box_leng);
    check_nearlist_correct(ptcls, near_ptcl_id, box_leng);

    // for debug
    // generate_graph_naive(ptcls, near_ptcl_id_naive, box_leng);
    // check(near_ptcl_id, near_ptcl_id_naive);
      
    reg_flag.resize(p_size, false);
    int   cur_patch = 0;
    search_patch(reg_flag, near_ptcl_id, patch_id, ptcls, id_stack, cur_patch);

    print_inout_elem(patch_id, fout, cur_patch, time, ptcls);

    if (time == 505000) {
      dump_one_frame(ptcls, patch_id);
      // dump_near_ptcl(ptcls, near_ptcl_id);
      calc_min_hei(ptcls, patch_id);
    }

    ptcl_id.clear();
    reg_flag.clear();
    
    for (size_t i = 0; i < near_ptcl_id.size(); i++)
      near_ptcl_id[i].clear();
    for (size_t i = 0; i < near_ptcl_id_naive.size(); i++)
      near_ptcl_id_naive[i].clear();
  }
}
