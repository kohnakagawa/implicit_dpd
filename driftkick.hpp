#pragma once

template<class Tpsys>
void KickAndDrift(Tpsys& sys,
		  const PS::F64 dt,
		  const PS::F64vec& box,
		  const PS::F64vec& ibox)
{
  PS::S64 n = sys.getNumberOfParticleLocal();
  for(int i = 0; i < n; i++) {
    sys[i].acc_bef = sys[i].acc;
    sys[i].vel_bef = sys[i].vel;
    
    const PS::F64 dr = dt * (sys[i].vel + 0.5 * dt * sys[i].acc) ;
    const PS::F64 r_temp = sys[i].pos + dr;

    sys[i].pos.x = r_temp.x - floor(r_temp.x * ibox.x) * box.x;
    sys[i].pos.y = r_temp.y - floor(r_temp.y * ibox.y) * box.y;
    sys[i].pos.z = r_temp.z - floor(r_temp.z * ibox.z) * box.z;

    //sys[i].vel += 0.5 * dt * sys[i].acc;
    sys[i].vel += 0.65 * dt * sys[i].acc;
  }
}

template<class Tpsys>
void Kick(Tpsys& sys,
	  const PS::F64 dt)
{
  PS::S64 n = sys.getNumberOfParticleLocal();
  for(int i = 0; i < n; i++) {
    sys[i].vel = sys[i].vel_bef + dt * 0.5 * (sys[i].acc + sys[i].acc_bef);
  }
}

