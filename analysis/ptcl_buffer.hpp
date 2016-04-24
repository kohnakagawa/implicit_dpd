#pragma once

#include <array>
#include "particle_simulator.hpp"
#include "io_util.hpp"
#include "parameter.hpp"
#include "ptcl_class.hpp"

struct Particle {
  PS::F64vec r;
  int prop, amp_id, unit;
  int hash = -1;

  // essential member functions
  void SetHash(const std::array<int, 3>& grid_numb,
	       const std::array<double, 3>& grid_leng) {
    int i = static_cast<int>(r[0] / grid_leng[0]), j = static_cast<int>(r[1] / grid_leng[1]), k = static_cast<int>(r[2] / grid_leng[2]);
 
    if (i == grid_numb[0]) i--;
    if (j == grid_numb[1]) j--;
    if (k == grid_numb[2]) k--;
  
    hash = i + grid_numb[0] * (j + k * grid_numb[1]);
    assert((hash >= 0 && hash < grid_numb[0] * grid_numb[1] * grid_numb[2]) && "hash value is out of range.");
  }
  
  int GetHash() const {
    return hash;
  }

  void CopyFromPtclOrg(const FPDPD& ptcl) {
    r = ptcl.pos;
    prop = ptcl.prop;
    amp_id = ptcl.amp_id;
    unit = ptcl.unit;
  }
  
  bool IsTarget() const {
    return ((prop == Parameter::Hyphil) || (prop == Parameter::Hyphob));
  }
};
