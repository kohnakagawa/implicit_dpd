#pragma once

#include "parameter.hpp"

template<class Tpsys>
void kick_and_drift(Tpsys& sys,
		    const PS::F64 dt,
		    const PS::F64vec& box,
		    const PS::F64vec& ibox)
{
  PS::S32 n = sys.getNumberOfParticleLocal();
  for(PS::S32 i = 0; i < n; i++) {
    const PS::F64vec dr = dt * (sys[i].vel + 0.5 * dt * sys[i].acc);
    const PS::F64vec r_temp = sys[i].pos + dr;

    sys[i].pos.x = r_temp.x - std::floor(r_temp.x * Parameter::ibox_leng.x) * Parameter::box_leng.x;
    sys[i].pos.y = r_temp.y - std::floor(r_temp.y * Parameter::ibox_leng.y) * Parameter::box_leng.y;
    sys[i].pos.z = r_temp.z - std::floor(r_temp.z * Parameter::ibox_leng.z) * Parameter::box_leng.z;

    sys[i].vel_buf = sys[i].vel + 0.5 * dt * sys[i].acc; //NOTE: vel_buf is needed for velocity update.

    //sys[i].vel += 0.5 * dt * sys[i].acc;
    sys[i].vel += 0.65 * dt * sys[i].acc;
  }
}

template<class Tpsys>
void kick(Tpsys& sys,
	  const PS::F64 dt)
{
  PS::S32 n = sys.getNumberOfParticleLocal();
  for(PS::S32 i = 0; i < n; i++)
    sys[i].vel = sys[i].vel_buf + dt * 0.5 * sys[i].acc;
}

