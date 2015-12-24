#pragma once

#include "parameter.hpp"

template<class Tpsys>
void drift_and_predict(Tpsys& sys,
		       const PS::F64 dt,
		       const PS::F64vec& box,
		       const PS::F64vec& ibox)
{
  PS::S32 n = sys.getNumberOfParticleLocal();
  for(PS::S32 i = 0; i < n; i++) {
    sys[i].pos += dt * (sys[i].vel + 0.5 * dt * sys[i].acc);

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

