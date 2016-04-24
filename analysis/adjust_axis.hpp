#pragma once

#include "Eigen/Dense"
#include "particle_simulator.hpp"

template<class Ptcl>
class AxisAdjuster {
  PS::F64vec cmpos;
  Eigen::Matrix3d Im;
  
  PS::F64vec CalcCmPos(const std::vector<Ptcl>& ptcls) {
    cmpos.x = cmpos.y = cmpos.z = 0.0;
    for (size_t i = 0; i < ptcls.size(); i++) {
      cmpos += ptcls[i].r;
    }
    cmpos /= ptcls.size();
    return cmpos;
  }
  
public:
  void CreateMomInertia(const std::vector<Ptcl>& ptcls) {
    CalcCmPos(ptcls);
    Im = Eigen::Matrix3d::Zero();
    for (size_t pi = 0; pi < ptcls.size(); pi++) {
      const auto cm2ri = ptcls[pi].r - cmpos;
      for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
	  Im(i, j) += cm2ri[i] * cm2ri[j];
    }

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
	Im(i, j) /= ptcls.size();
  }

  void DoTransform(std::vector<Ptcl>& ptcls) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(Im);
    if (es.info() != Eigen::Success) {
      std::cerr << "Error occurs in Eigen solver.\n";
      std::exit(1);
    }
    
    Eigen::Vector3d::Index maxId;
    const auto max = es.eigenvalues().maxCoeff(&maxId);
    const auto evec = es.eigenvectors();
    
    Eigen::Vector3d new_base[3];
    int id = 0;
    for (int i = 0; i < 3; i++) {
      if (i != maxId) {
	new_base[id++] = evec[i];
      } else {
	new_base[2] = evec[maxId];
      }
    }
    
    // change base
    Eigen::Vector3d cm2ptcl;
    for (size_t pi = 0; pi < ptcls.size(); pi++) {
      cm2ptcl << (ptcls[pi].r[0] - cmpos[0]), (ptcls[pi].r[1] - cmpos[1]), (ptcls[pi].r[2] - cmpos[2]);
      ptcls[pi].r.x = cm2ptcl * new_base[0];
      ptcls[pi].r.y = cm2ptcl * new_base[1];
      ptcls[pi].r.z = cm2ptcl * new_base[2];
    }
  }
};
