#pragma once
#include "trjanalysis.hpp"

class TrjAnalysisThickness : public TrjAnalysis<FPDPD, Particle> {
  PS::S32 GetGridIdx(const PS::F64 rx,
		     const PS::F64 ry,
		     const PS::F64 grid_len_x,
		     const PS::F64 grid_len_y,
		     const std::vector<PS::S32>& grid) {
    PS::S32 idx[] = {
      static_cast<PS::S32>(rx / grid_len_x),
      static_cast<PS::S32>(ry / grid_len_y)
    };
    if (idx[0] == grid[0]) idx[0]--;
    if (idx[1] == grid[1]) idx[1]--;
    
    return idx[0] + idx[1] * grid[0];
  }

  PS::F64 CalcLocalThickness(std::vector<std::pair<PS::S32, PS::S32> >& cnt,
			     std::vector<std::pair<PS::F64, PS::F64> >& height_in_grid,
			     std::vector<std::pair<PS::F64, PS::F64> >& protru_in_grid,
			     std::vector<PS::F64>& thick_in_grid,
			     const std::vector<PS::S32>& patch_id,
			     const std::vector<PS::S32>& grid,
			     const std::vector<PS::F64>& grid_len,
			     const std::vector<PS::S32>& laxis,
			     const PS::S32 naxis) {
    const auto num_ptcls = ptcls_.size();
    for (auto i = 0u; i < num_ptcls; i++) {
      const auto hash = GetGridIdx(ptcls_[i].r[laxis[0]],
				   ptcls_[i].r[laxis[1]],
				   grid_len[laxis[0]],
				   grid_len[laxis[1]],
				   grid);
      if (patch_id[i] == 0) {
	height_in_grid[hash].first += ptcls_[i].r[naxis];
	cnt[hash].first++;
      } else {
	height_in_grid[hash].second += ptcls_[i].r[naxis];
	cnt[hash].second++;
      }
    }
    
    for (size_t i = 0; i < height_in_grid.size(); i++) {
      height_in_grid[i].first  /= cnt[i].first;
      height_in_grid[i].second /= cnt[i].second;
      thick_in_grid[i] = std::abs(height_in_grid[i].first - height_in_grid[i].second);
    }
    
    for (auto i = 0u; i < num_ptcls; i++) {
      const auto hash = GetGridIdx(ptcls_[i].r[laxis[0]],
				   ptcls_[i].r[laxis[1]],
				   grid_len[laxis[0]],
				   grid_len[laxis[1]],
				   grid);
      if (patch_id[i] == 0) {
	const auto dist = (ptcls_[i].r[naxis] - height_in_grid[hash].first);
	protru_in_grid[hash].first += dist * dist;
      } else {
	const auto dist = (ptcls_[i].r[naxis] - height_in_grid[hash].second);
	protru_in_grid[hash].second += dist * dist;
      }
    }

    for (size_t i = 0; i < protru_in_grid.size(); i++) {
      protru_in_grid[i].first  /= cnt[i].first;
      protru_in_grid[i].second /= cnt[i].second;
      protru_in_grid[i].first  = std::sqrt(protru_in_grid[i].first);
      protru_in_grid[i].second = std::sqrt(protru_in_grid[i].second);
    }
  }
  
public:
  TrjAnalysisThickness(const std::string cur_dir,
		 const char* trj_fname,
		 const char* param_fname) : TrjAnalysis(cur_dir, trj_fname, param_fname) {}
  
  ~TrjAnalysisThickness() override {}
  
  void DoAnalysis() override {
    const int naxis = 1;
    std::cerr << "We assume normal axis is y.\n";

    SetSearchRadius(1.1, 1.05);
    ptr_connector->Initialize(est_grid_leng_, cutof_leng_, Parameter::box_leng);
    
    constexpr PS::F64 init_grid_len = 3.0;

    std::vector<PS::S32> grid;
    std::vector<PS::F64> grid_len = {init_grid_len, init_grid_len};
    std::vector<PS::S32> laxis;
    for (int i = 0; i < 3; i++) {
      if (i != naxis) {
	const auto num = static_cast<int>(Parameter::box_leng[i] / init_grid_len);
	grid.push_back(num);
	grid_len.push_back(Parameter::box_leng[i] / num);
	laxis.push_back(i);
      }
    }
    
    const auto all_grid = grid[0] * grid[1];
    std::vector<std::pair<PS::S32, PS::S32> > cnt_in_grid(all_grid, {0, 0});
    std::vector<std::pair<PS::F64, PS::F64> > height_in_grid(all_grid, {0.0, 0.0});
    std::vector<std::pair<PS::F64, PS::F64> > protru_in_grid(all_grid, {0.0, 0.0});
    std::vector<PS::F64> thick_in_grid(all_grid, 0.0);
    
    std::string fname = cur_dir_ + "/thickness_protoru.txt";
    std::ofstream fout(fname.c_str());
    fout << std::setprecision(15);
   
    PS::U32 time = 0;
    while (true) {
      if (ReadOneFrame(time, [] (const Particle& p) -> bool {return p.prop == Parameter::Hyphil;})) break;
      std::cout << "time = " << time << std::endl;
      ptr_connector->ResizeIfNeeded(ptcls_.size());

      const auto num_patch = gen_connected_image(ptr_connector.get(), ptcls_, Parameter::box_leng, true);
      std::cout << "# of patch is " << num_patch << std::endl;

      if (num_patch == 2) {
	CalcLocalThickness(cnt_in_grid,
			   height_in_grid,
			   protru_in_grid,
			   thick_in_grid,
			   ptr_connector->patch_id(),
			   grid,
			   grid_len,
			   laxis,
			   naxis);
	fout << std::accumulate(protru_in_grid.cbegin(),
				protru_in_grid.cend(),
				0.0,
				[] (const double sum, const std::pair<PS::F64, PS::F64>& elem) {
				    return sum + elem.first + elem.second;
				  }) / all_grid
	     << " ";
	fout << std::accumulate(thick_in_grid.cbegin(), thick_in_grid.cend(), 0.0) / all_grid << std::endl;
	cnt_in_grid.assign(cnt_in_grid.size(), {0, 0});
	height_in_grid.assign(height_in_grid.size(), {0, 0});
	protru_in_grid.assign(height_in_grid.size(), {0.0, 0.0});
	thick_in_grid.assign(thick_in_grid.size(), 0.0);
      }
      
      ptr_connector->ClearForNextStep();      
    }
  }
};
