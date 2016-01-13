#pragma once

struct FPDPD {
  PS::U32 id, id_loc, prop, amp_id, unit;
  PS::F64vec pos, delta_sumr;
  PS::F64vec vel, vel_buf;
  PS::F64vec acc;
  PS::F64vec press;
  std::array<PS::F64, Parameter::prop_num> density;

  //essential member functions
  PS::F64 getRSearch() const {
    return Parameter::search_rad;
  }
  PS::F64vec getPos() const {
    return this->pos;
  }
  void setPos(const PS::F64vec& p) {
    pos = p;
  }
  
#ifdef USE_NEIGH_BUFFER
  void copyFromForce(const RESULT::ForceDPDDensity& force) {
    acc		= force.acc;
    press	= force.press;
    density     = force.dens;
  }  
#else
  void copyFromForce(const RESULT::ForceDPD& force) {
    acc		= force.acc;
    press	= force.press;
  }
  void copyFromForce(const RESULT::Density& dens) {
    density = dens.dens;
  }
#endif
  
  //for I/O
  void readAscii(FILE *fp) {
    char buf;
    fscanf(fp, "%c %lf %lf %lf %u %u %u %u %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	   &buf, &(pos.x), &(pos.y), &(pos.z),
	   &id, &prop, &amp_id, &unit,
	   &(vel.x), &(vel.y), &(vel.z), &(vel_buf.x), &(vel_buf.y), &(vel_buf.z),
	   &(acc.x), &(acc.y), &(acc.z));
    delta_sumr = 0.0;
  }
  void writeAscii(FILE *fp) const {
    fprintf(fp, "%c %.15g %.15g %.15g %u %u %u %u %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n",
	    Parameter::atom_type[prop], pos.x, pos.y, pos.z,
	    id, prop, amp_id, unit,
	    vel.x, vel.y, vel.z, vel_buf.x, vel_buf.y, vel_buf.z,
	    acc.x, acc.y, acc.z);
  }
};
