#pragma once

#include "parameter.hpp"

template<class Tpsys>
void drift_and_predict(Tpsys& sys,
		       const PS::F64 dt,
		       const PS::F64vec& box_len,
		       const PS::F64vec& ibox_len)
{
  PS::S32 n = sys.getNumberOfParticleLocal();
  for(PS::S32 i = 0; i < n; i++) {
    const PS::F64vec dr(dt * (sys[i].vel.x + 0.5 * dt * sys[i].acc.x),
			dt * (sys[i].vel.y + 0.5 * dt * sys[i].acc.y),
			dt * (sys[i].vel.z + 0.5 * dt * sys[i].acc.z) );
    
    const PS::F64vec temp_r = sys[i].pos + dr;

    sys[i].delta_sumr += dr;

    sys[i].pos.x = temp_r.x - std::floor(temp_r.x * ibox_len.x) * box_len.x;
    sys[i].pos.y = temp_r.y - std::floor(temp_r.y * ibox_len.y) * box_len.y;
    sys[i].pos.z = temp_r.z - std::floor(temp_r.z * ibox_len.z) * box_len.z;

    sys[i].vel_buf = sys[i].vel + 0.5 * dt * sys[i].acc; //NOTE: vel_buf is needed for velocity update.

    sys[i].vel += 0.5 * dt * sys[i].acc; //NOTE: lambda = 0.5 is better than lambda = 0.65 when dt = 0.005.
    //sys[i].vel += 0.65 * dt * sys[i].acc;
  }
}

template<class Tpsys>
void kick(Tpsys& sys,
	  const PS::F64 dt)
{
  PS::S32 n = sys.getNumberOfParticleLocal();
  for (PS::S32 i = 0; i < n; i++) {
    sys[i].vel = sys[i].vel_buf + dt * 0.5 * sys[i].acc;
  }
}

template <class Tpsys>
void remove_cmdrift_global(Tpsys& sys) {
  PS::F64vec vel_sum_loc(0.0, 0.0, 0.0);
  const PS::S32 n_loc = sys.getNumberOfParticleLocal();
  for (PS::S32 i = 0; i < n_loc; i++) {
    vel_sum_loc += sys[i].vel;
  }
  const PS::S32 n = sys.getNumberOfParticleGlobal();
  const PS::F64vec vel_sum = PS::Comm::getSum(vel_sum_loc) / PS::F64(n);
  for (PS::S32 i = 0; i < n_loc; i++) {
    sys[i].vel -= vel_sum;
  }
}
