#pragma once

#include "Eigen/Eigenvalues"
#include "particle_simulator.hpp"

template<class Ptcl>
class AxisAdjuster {
  PS::F64vec cmpos;
  Eigen::Matrix3d Im;
  
  PS::F64vec CalcCmPos(const std::vector<Ptcl>& ptcls,
		       const std::vector<int>& patch_id,
		       const int tar_ptch_id) {
    int cnt = 0;
    cmpos.x = cmpos.y = cmpos.z = 0.0;
    for (size_t i = 0; i < ptcls.size(); i++) {
      if (patch_id[i] == tar_ptch_id) {
	cmpos += ptcls[i].r;	
	cnt++;
      }
    }
    cmpos /= cnt;
    return cmpos;
  }
  
public:
  void CreateMomInertia(const std::vector<Ptcl>& ptcls,
			const std::vector<int>& patch_id,
			const int tar_ptch_id) {
    CalcCmPos(ptcls, patch_id, tar_ptch_id);
    Im = Eigen::Matrix3d::Zero();
    int cnt = 0;
    for (size_t pi = 0; pi < ptcls.size(); pi++) {
      if (patch_id[pi] == tar_ptch_id) {
	const auto cm2ri = ptcls[pi].r - cmpos;
	for (int i = 0; i < 3; i++)
	  for (int j = 0; j < 3; j++)
	    Im(i, j) += cm2ri[i] * cm2ri[j];
	cnt++;
      }
    }

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
	Im(i, j) /= cnt;
  }

  void DoTransform(std::vector<Ptcl>& ptcls) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(Im);
    if (es.info() != Eigen::Success) {
      std::cerr << "Error occurs in Eigen solver.\n";
      std::exit(1);
    }
    
    Eigen::Vector3d::Index maxId;
    const auto max = es.eigenvalues().maxCoeff(&maxId);
    std::cout << "max eigen value is " << max << std::endl;
    const auto evec = es.eigenvectors();
    
    Eigen::Vector3d new_base[3];
    int id = 0;
    for (int i = 0; i < 3; i++) {
      if (i != maxId) {
	new_base[id++] = evec.col(i);
      } else {
	new_base[2] = evec.col(maxId);
      }
    }
    
    // change base
    Eigen::Vector3d cm2ptcl;
    for (size_t pi = 0; pi < ptcls.size(); pi++) {
      cm2ptcl << (ptcls[pi].r[0] - cmpos[0]), (ptcls[pi].r[1] - cmpos[1]), (ptcls[pi].r[2] - cmpos[2]);
      ptcls[pi].r.x = cm2ptcl.dot(new_base[0]);
      ptcls[pi].r.y = cm2ptcl.dot(new_base[1]);
      ptcls[pi].r.z = cm2ptcl.dot(new_base[2]);
    }
  }
};
