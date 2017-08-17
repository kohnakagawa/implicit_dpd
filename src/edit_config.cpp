#define _FILE_OFFSET_BITS 64 // This line is needed on 32bit environment.

#include <iostream>
#include <string>
#include <vector>
#include <cstdio>
#include <limits>
#include "particle_simulator.hpp"
#include "io_util.hpp"
#include "parameter.hpp"
#include "f_calculator.hpp"

PS::F64vec Parameter::box_leng, Parameter::ibox_leng;
PS::U32 Parameter::all_time, Parameter::step_mic, Parameter::step_mac, Parameter::time;

constexpr char Parameter::atom_type[21];
const char axis_name[] = {
  'X', 'Y', 'Z'
};

PS::F64 gen_thermal_vel(const PS::F64 Tempera) {
  const PS::F64 uni_rand[] = {
    PS::MT::genrand_res53(),
    PS::MT::genrand_res53()};
  return std::sqrt(-2.0 * Tempera * std::log(uni_rand[0])) * std::cos(2.0 * M_PI * uni_rand[1]);
}

PS::F64vec gen_thermal_vel() {
  const PS::F64 T = 1.0;
  return PS::F64vec(gen_thermal_vel(T),
		    gen_thermal_vel(T),
		    gen_thermal_vel(T));
}

#define CHECK_DIM(dim, axis)                                  \
  do {                                                        \
    if (dim == 2) {                                           \
      assert(axis < 3);                                       \
    } else if (dim == 3) {                                    \
      assert(axis == 0xffffffff);                             \
    } else {                                                  \
      std::cerr << "Please check dimenstion.\n";              \
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;	\
      std::exit(1);                                           \
    }                                                         \
  } while(false)

PS::F64vec gen_pos(const PS::F64vec& leng) {
  return PS::F64vec(leng.x * PS::MT::genrand_res53(),
		    leng.y * PS::MT::genrand_res53(),
		    leng.z * PS::MT::genrand_res53());
}

void append_solvent(std::vector<FPDPD>& ptcls,
		    const PS::U32 new_id,
		    const PS::F64vec& pos,
		    const PS::F64vec& vel) {
    FPDPD solvent;
    solvent.id = new_id;
    solvent.prop = Parameter::Solvent;
    solvent.amp_id = solvent.unit = 0xffffffff;
    solvent.pos = pos;
    solvent.vel = solvent.vel_buf = vel;
    solvent.acc = solvent.press = 0.0;
    ptcls.push_back(solvent);
}

template<class FP>
void add_solvents_cuboid(std::vector<FP>& ptcls,
			 const PS::F64vec& up_pos,
			 const PS::F64vec& dw_pos,
			 const PS::U32 num) {
  const PS::F64vec leng = up_pos - dw_pos;
  PS::U32 added_num = 0, new_id = ptcls.size();
  while (added_num < num) {
    const PS::F64vec pos = gen_pos(leng) + dw_pos;
    const PS::F64vec vel = gen_thermal_vel();
    append_solvent(ptcls, new_id++, pos, vel);
    added_num++;
  }
}

template<class FP>
void add_solvents_sph(std::vector<FP>& ptcls,
		      const PS::F64vec& cent,
		      const PS::F64 rad,
		      const PS::U32 num,
		      const PS::F64vec& box_leng,
		      const bool add_in,
		      const PS::U32 dim,
		      const PS::U32 axis) {
  CHECK_DIM(dim, axis);

  PS::U32 added_num = 0, new_id = ptcls.size();
  while (added_num < num) {
    const PS::F64vec pos = gen_pos(box_leng);
    PS::F64vec cent2pos_vec = pos - cent;
    PS::F64 cent2pos2 = cent2pos_vec * cent2pos_vec;
    if (dim == 2) {
      cent2pos2 -= cent2pos_vec[axis] * cent2pos_vec[axis];
    }
    const auto cent2pos = std::sqrt(cent2pos2);
    const bool is_in = (cent2pos < rad);
    if (is_in ^ (!add_in)) {
      const PS::F64vec vel = gen_thermal_vel();
      append_solvent(ptcls, new_id++, pos, vel);
      added_num++;
    }
  }
}

template<class FP>
void remove_cm_drift(std::vector<FP>& ptcls) {
  PS::F64vec cm_vel = 0.0;
  for (const auto& ptcl : ptcls) {
    cm_vel += ptcl.vel;
  }
  cm_vel /= ptcls.size();

  for (auto& ptcl : ptcls) {
    ptcl.vel -= cm_vel;
  }
}

// Func = std::max_element or std::min_element
template<class FP, class Func>
PS::U32 guess_bilayer_axis(const std::vector<FP>& ptcls,
			   PS::F64& memb_dw,
			   PS::F64& memb_up,
			   Func func) {
  PS::F64vec downer = std::numeric_limits<PS::F64>::max(), upper = 0.0;
  for (const auto& ptcl : ptcls) {
    upper  = ApplyEach([](const PS::F64& v0, const PS::F64& v1) { return std::max(v0, v1); }, ptcl.pos, upper);
    downer = ApplyEach([](const PS::F64& v0, const PS::F64& v1) { return std::min(v0, v1); }, ptcl.pos, downer);
  }
  const PS::F64vec leng = upper - downer;
  const auto ptr_beg = &(leng[0]);
  const PS::U32 normal_axis = std::distance(ptr_beg, func(ptr_beg, ptr_beg + 3));
  memb_up = upper[normal_axis];
  memb_dw = downer[normal_axis];
  return normal_axis;
}

template<class FP>
void guess_bilayer_cent_and_rad(const std::vector<FP>& ptcls,
				PS::F64& rad,
				PS::F64vec& cent,
				const PS::U32 dim,
				const PS::U32 axis) {
  // ASSUME: spherical vesicle is located in box center.
  // NOTE: Do not check a connectivity of this vesicle.
  CHECK_DIM(dim, axis);

  cent = 0.0;
  for (const auto& ptcl : ptcls) {
    cent += ptcl.pos;
  }
  cent /= ptcls.size();

  rad = 0.0;
  for (const auto& ptcl : ptcls) {
    const auto cent2ptcl = ptcl.pos - cent;
    auto rad2 = cent2ptcl * cent2ptcl;
    if (dim == 2) rad2 -= cent2ptcl[axis] * cent2ptcl[axis];
    rad += std::sqrt(rad2);
  }
  rad /= ptcls.size();
}

template<class FP>
void check_config_is_valid(std::vector<FP>& ptcls, const PS::F64vec& box) {
  PS::F64 Tempera = 0.0;
  for (const auto& ptcl : ptcls) {
    for (int i = 0; i < 3; i++)
      assert(ptcl.pos[i] >= 0.0 && ptcl.pos[i] <= box[i]);
    if (ptcl.prop != Parameter::Solvent) {
      assert(ptcl.amp_id != 0xffffffff);
      assert(ptcl.unit != 0xffffffff);
    } else {
      assert(ptcl.amp_id == 0xffffffff);
      assert(ptcl.unit == 0xffffffff);
    }

    Tempera += (ptcl.vel * ptcl.vel) / 3.0;
  }
  Tempera /= ptcls.size();
  std::cerr << "System temperature is " << Tempera << std::endl;
}

template<class FP>
bool read_one_frame(std::vector<FP>& ptcls,
		    PS::U32& num,
		    PS::U32& time,
		    FILE*& fp) {
  char name[4];
  assert(0 != fscanf(fp, "%u\n", &num));
  if (feof(fp)) return false;
  assert(0 != fscanf(fp, "%s %u\n", name, &time));
  ptcls.resize(num);
  for (auto& ptcl : ptcls) ptcl.readAscii(fp);
  return true;
}

template<class FP>
void read_config(std::vector<FP>& ptcls,
		 PS::U32& num,
		 PS::U32& time,
		 const char* fname) {
  FILE* fp = io_util::xfopen(fname, "r");
  assert(read_one_frame(ptcls, num, time, fp));
  char buf;
  assert(0 != fscanf(fp, "%c\n", &buf));
  if (!feof(fp)) {
    std::cerr << "# of lines is not equal to the information specified in file header.\n";
    std::exit(1);
  }
  fclose(fp);
}

template<class FP>
void read_trj(std::vector<FP>& ptcls,
	      PS::U32& num,
	      PS::U32& time,
	      const PS::U32 tar_time,
	      const char* fname) {
  FILE* fp = io_util::xfopen(fname, "r");
  while (true) {
    if (!read_one_frame(ptcls, num, time, fp)) break;
    std::cout << "reading frame at time = " << time << std::endl;
    if (tar_time == time) {
      std::cout << "Found!.\n";
      break;
    }
  }
  fclose(fp);

  if (time != tar_time) {
    std::cerr << "Error!\n";
    std::cerr << "Tried to find the frame at time = " << tar_time << ".\n";
    std::cerr << "However, we can not find it.\n";
    std::exit(1);
  }
}

template<class FP>
void search_near_ptcl_id(std::vector<PS::U32>& near_ids,
			 std::vector<FP>& near_ptcls,
			 const std::vector<FP>& ptcls,
			 const PS::F64 rad,
			 const PS::S32 core_ptcl_id) {
  near_ids.reserve(100);
  near_ptcls.reserve(100);

  const PS::F64vec tar_pos = ptcls[core_ptcl_id].pos;
  const PS::F64 rad2 = rad * rad;

  const PS::S32 ptcl_size = ptcls.size();

  for (PS::S32 i = 0; i < ptcl_size; i++) {
    const PS::S32 prop = ptcls[i].prop;
    PS::F64vec core2ptcl = ptcls[i].pos - tar_pos;
    ForceBonded<PS::ParticleSystem<FPDPD> >::MinImage(core2ptcl);
    const PS::F64 core2ptcl_dist2 = core2ptcl * core2ptcl;

    if ((core2ptcl_dist2 < rad2) && (prop == Parameter::Hyphil)) {
      near_ids.push_back(ptcls[i].amp_id);
      near_ptcls.push_back(ptcls[i]);
    }
  }
}

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cerr << "argv[1] is xyz file name.\n";
    std::cerr << "argv[2] is target directory name.\n";
    std::exit(1);
  }
  const std::string fname = argv[1];
  const std::string cur_dir = argv[2];
  PS::MT::init_genrand(100);

  Parameter param(cur_dir);
  param.Initialize();
  param.LoadParam();

  std::string mode;
  std::cout << "Choose execution mode [solvent/search/trjextract/extractcore].\n";
  std::cin >> mode;

  std::vector<FPDPD> ptcls; PS::U32 num = 0, time = 0;
  if (mode != "trjextract") {
    // read only one frame
    read_config(ptcls, num, time, fname.c_str());
    std::cout << "Successfully load particle data.\n";
  }

  if (mode == "solvent") {
    std::string shape;
    std::cout << "Choose membrane shape [flat/sphere/cylind].\n";
    std::cin >> shape;
    if ((shape != "flat") && (shape != "sphere") && (shape != "cylind")) {
      std::cerr << "Choose correct membrane shape.\n";
      std::exit(1);
    }

    int added_sol_num = -1;
    std::cout << "Add solvent particles in " << fname << std::endl;
    std::cout << "# of solvent particles will be added.\n";
    std::cin  >> added_sol_num;

    if (added_sol_num <= 0) {
      std::cerr << "added_sol_num should be larger than 0.\n";
      std::exit(1);
    }

    const PS::F64 est_memb_thick = Parameter::all_unit * Parameter::bond_leng;

    if (shape == "flat") {
      std::cerr << "Choose upper side / downer side [up/dw].\n";
      std::string side;
      std::cin >> side;
      const bool b_up = (side == "up");

      PS::F64 memb_up = 0.0, memb_dw = 0.0;
      const auto axis = guess_bilayer_axis(ptcls, memb_dw, memb_up,
					   [](const PS::F64* beg, const PS::F64* end) { return std::min_element(beg, end); });
      std::cerr << "Bilayer normal axis is " << axis_name[axis] << ".\n";

      PS::F64vec up_pos = Parameter::box_leng, dw_pos = 0.0;
      if (b_up) {
        dw_pos[axis] = memb_up;
      } else {
        up_pos[axis] = memb_dw;
      }
      add_solvents_cuboid(ptcls, up_pos, dw_pos, added_sol_num);
    } else if ((shape == "sphere") || (shape == "cylind")) {
      PS::U32 dim = 3, axis = 0xffffffff;
      if (shape == "cylind") {
        dim = 2;
        PS::F64 memb_up = 0.0, memb_dw = 0.0;
        axis = guess_bilayer_axis(ptcls, memb_dw, memb_up,
                                  [](const PS::F64* beg, const PS::F64* end) { return std::max_element(beg, end);	});
      }
      PS::F64 rad = 0.0;
      PS::F64vec cent = 0.0;
      guess_bilayer_cent_and_rad(ptcls, rad, cent, dim, axis);

      std::cerr << "Vesicle radius is " << rad << std::endl;
      std::cerr << "Center position is " << cent << std::endl;

      std::string side;
      std::cerr << "Choose inside or outside [in/out].\n";
      std::cin >> side;
      const bool b_in = (side == "in");

      if (b_in) {
        rad -= est_memb_thick;
      } else {
        rad += est_memb_thick;
      }

      add_solvents_sph(ptcls, cent, rad, added_sol_num, Parameter::box_leng, b_in, dim, axis);
    }

    remove_cm_drift(ptcls);
    check_config_is_valid(ptcls, Parameter::box_leng);

    const std::string out_fname = cur_dir + "/rev_config.xyz";
    FILE* fp = io_util::xfopen(out_fname.c_str(), "w");
    io_util::WriteXYZForm(&(ptcls[0]), ptcls.size(), time, fp);
    fclose(fp);
    fp = nullptr;

    std::cout << "revised configuration is generaged at " << out_fname << std::endl;
  } else if (mode == "search") {
    PS::S32 core_ptcl_id = -1;
    std::cout << "Choose core ptcl id\n";
    std::cin >> core_ptcl_id;
    CHECK_GE(core_ptcl_id, 0);
    CHECK_LT(core_ptcl_id, static_cast<int>(ptcls.size()));
    CHECK_EQ(core_ptcl_id % Parameter::all_unit, 0);

    PS::F64 s_rad = 0.0;
    std::cout << "Now search particle ids which are located near id " << core_ptcl_id << ".\n";
    std::cout << "Search radius ?\n";
    std::cin >> s_rad;
    CHECK_GT(s_rad, 0.0);

    std::vector<PS::U32> near_ids;
    std::vector<FPDPD> near_ptcls;
    search_near_ptcl_id(near_ids, near_ptcls, ptcls, s_rad, core_ptcl_id);

    std::cout << "Detected # of near particle ids is " << near_ids.size() << std::endl << std::endl;;
    const std::string fname_conf = cur_dir + "/ptcl_ids_near_core.xyz";
    const std::string fname_ids  = cur_dir + "/nearid_lists.txt";
    FILE* fp = io_util::xfopen(fname_conf.c_str(), "w");
    io_util::WriteXYZForm(&(near_ptcls[0]), near_ptcls.size(), 0, fp);
    fclose(fp);
    fp = nullptr;

    std::ofstream fout(fname_ids.c_str());
    for (PS::U32 i = 0; i < near_ids.size(); i++) fout << near_ids[i] << " ";
    fout << std::endl;

    std::cout << "The configuration of these particles are written in " << fname_conf << std::endl;
    std::cout << "Particle ids near core amphiphile are written in " << fname_ids << std::endl << std::endl;;
  } else if (mode == "trjextract") {
    std::cout << "Please specify simulation time at which configuration data is extracted.\n";
    PS::U32 tar_time = 0;
    std::cin >> tar_time;
    std::cout << "Now read trajectory data.\n";
    read_trj(ptcls, num, time, tar_time, fname.c_str());

    std::stringstream ss;
    ss << cur_dir << "/extracted_frame_time" << tar_time << ".xyz";
    FILE* fp = io_util::xfopen(ss.str().c_str(), "w");
    io_util::WriteXYZForm(&(ptcls[0]), num, tar_time, fp);
    std::cout << "Extracted data -> " << ss.str() << std::endl;
    fclose(fp);
    fp = nullptr;

  } else if (mode == "extractcore") {
    const std::string fname = cur_dir + "/core_amp_coord.xyz";

    std::vector<FPDPD> core_ptcls;
    for (const auto& ptcl : ptcls)
      for (auto j = 0u; j < param.core_amp_id_.size(); j++)
        if (ptcl.amp_id == param.core_amp_id_[j])
          core_ptcls.push_back(ptcl);

    FILE* fp = io_util::xfopen(fname.c_str(), "w");
    io_util::WriteXYZForm(&(core_ptcls[0]), core_ptcls.size(), 0, fp);
    fclose(fp);

    std::cout << "The configuration of these particles are written in " << fname << "\n";
  } else {
    std::cerr << "I do not know what to do.\n";
    std::exit(1);
  }
}
