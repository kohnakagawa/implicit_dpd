#include <iostream>
#include <algorithm>
#include "particle_simulator.hpp"
#include "io_util.hpp"
#include "parameter.hpp"
#include "f_calculator.hpp"

const PS::F64vec box(30.0);

struct less {
  bool operator() (const FPDPD& p1, const FPDPD& p2) const {
    const int key1 = p1.amp_id * Parameter::all_unit + p1.unit;
    const int key2 = p2.amp_id * Parameter::all_unit + p2.unit;
    return key1 < key2;
  }
};

void MinImage(PS::F64vec& drij) {
  drij.x -= box.x * std::round(drij.x / box.x);
  drij.y -= box.y * std::round(drij.y / box.y);
  drij.z -= box.z * std::round(drij.z / box.z);
}

int main(int argc, char *argv[]) {
  if(argc != 3) {
    std::cerr << "argv[1] is target directory.\n";
    std::cerr << "argv[2] is beg time\n";
    PS::Abort();
  }
  
  std::vector<FPDPD> ptcls;
  const std::string fname = std::string(argv[1]) + "/traject.xyz";
  const int beg_time = std::atoi(argv[2]);
  
  FILE* fp = fopen(fname.c_str(), "r");
  if(fp == nullptr) std::exit(1);

  double mean_bond_all = 0.0;
  double mean_angle_all = 0.0;
  int cnt = 0;

  FILE* fout_bht = fopen("bondht_hist.txt", "w");
  FILE* fout_btt = fopen("bondtt_hist.txt", "w");
  FILE* fout_a = fopen("angle_hist.txt", "w");
  FILE* fout_e = fopen("h2e.txt", "w");
  
  while(true) {
    int cur_time, num;
    char comment[4] = {'\0'};
    fscanf(fp, "%d\n", &num);
    fscanf(fp, "%s %d\n", comment, &cur_time);
  
    std::cout << comment << " " << cur_time << std::endl;

    
    if(feof(fp)) break;
    
    ptcls.resize(num);
    for(int i = 0; i < num; i++)
      ptcls[i].readAscii(fp);
    
    std::sort(ptcls.begin(), ptcls.end(), less());

    double mean_bond = 0.0, mean_angle = 0.0;
    for(int pi = 0; pi < num; pi += Parameter::all_unit) {
      PS::F64vec h2e_vec = ptcls[pi + Parameter::all_unit - 1].pos - ptcls[pi].pos;
      MinImage(h2e_vec);
      const double h2e_vec_norm = std::sqrt(h2e_vec * h2e_vec);
      if(cur_time >= beg_time)
	fprintf(fout_e, "%.15g\n", h2e_vec_norm);
      
      PS::F64vec ini_vec = ptcls[pi + 1].pos - ptcls[pi].pos;
      MinImage(ini_vec);
      mean_bond = std::sqrt(ini_vec * ini_vec);
      if(cur_time >= beg_time)
	fprintf(fout_bht, "%.15g\n", mean_bond);
      mean_angle = 0.0;
    
      for(PS::U32 i = 2; i < Parameter::all_unit; i++) {
      	PS::F64vec dr0 = ptcls[pi + i - 1].pos - ptcls[pi + i - 2].pos;
      	PS::F64vec dr1 = ptcls[pi + i    ].pos - ptcls[pi + i - 1].pos;

      	MinImage(dr0); MinImage(dr1);

      	const double dr0_norm = std::sqrt(dr0 * dr0);
      	const double dr1_norm = std::sqrt(dr1 * dr1);
      
      	double cos = dr0 * dr1 / (dr0_norm * dr1_norm );
      	if(cos > 1.0) cos = 1.0;
      	if(cos < -1.0) cos = -1.0;
      	const double theta = std::acos(cos);

      	assert(theta >= 0.0 && theta <= M_PI);
      
      	mean_angle += theta;
      	mean_bond  += dr1_norm;
	if(cur_time >= beg_time){
	  fprintf(fout_btt, "%.15g\n", dr1_norm);
	  fprintf(fout_a, "%.15g\n", theta);
	}
      }
      mean_bond /= (Parameter::all_unit - 1);
      mean_angle /= (Parameter::all_unit - 2);
    }
    
    if(cur_time >= beg_time) {
      mean_bond_all += mean_bond;
      mean_angle_all += mean_angle;
      cnt++;
    }
  }
  mean_angle_all /= cnt;
  mean_bond_all /= cnt;
  
  std::cout << "mean_bond_all" << std::endl;
  std::cout << mean_bond_all << std::endl;
  std::cout << "mean_angle_all" << std::endl;
  std::cout << mean_angle_all << std::endl;
  
  fclose(fp);
  fclose(fout_bht);
  fclose(fout_btt);
  fclose(fout_a);
  fclose(fout_e);
}
