#pragma once

#include <array>

namespace RESULT {
  struct ForceDPD {
    PS::F64vec acc;
    PS::F64vec press;
  
    void clear(){
      acc   = 0.0;
      press = 0.0;
    }
  };

  struct Density {
    std::array<PS::F64, Parameter::prop_num> dens;
    
    void clear() {
      dens.fill(0.0);
    }
  };
};

struct FPDPD {
  PS::U32 id, prop, amp_id, unit;
  PS::F64vec pos;
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
  void copyFromForce(const RESULT::ForceDPD& force) {
    acc = force.acc;
    press = force.press;
  }

  void copyFromForce(const RESULT::Density& dens) {
    density = dens.dens;
  }
  
  //for I/O
  void readAscii(FILE *fp) {
    char buf;
    fscanf(fp, "%c %lf %lf %lf %u %u %u %u %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	   &buf, &(pos.x), &(pos.y), &(pos.z),
	   &id, &prop, &amp_id, &unit,
	   &(vel.x), &(vel.y), &(vel.z), &(vel_buf.x), &(vel_buf.y), &(vel_buf.z),
	   &(acc.x), &(acc.y), &(acc.z));

  }
  void writeAscii(FILE *fp) const {
    fprintf(fp, "%c %.15g %.15g %.15g %u %u %u %u %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n",
	    Parameter::atom_type[prop], pos.x, pos.y, pos.z,
	    id, prop, amp_id, unit,
	    vel.x, vel.y, vel.z, vel_buf.x, vel_buf.y, vel_buf.z,
	    acc.x, acc.y, acc.z);
  }
};

namespace EPI {
  struct DPD {
    PS::U32 id, prop;
    PS::F64vec pos, vel;
    std::array<PS::F64, Parameter::prop_num> dens;

    PS::F64vec getPos() const {
      return this->pos;
    }

    PS::F64 getRSearch() const {
      return Parameter::search_rad;
    }

    void copyFromFP(const FPDPD& fp) {
      this->pos		= fp.pos;
      this->vel		= fp.vel;
      this->id		= fp.id;
      this->prop	= fp.prop;
      this->dens        = fp.density;
    }

  };
  
  struct Density {
    PS::U32 id, prop;
    PS::F64vec pos;

    PS::F64vec getPos() const {
      return this->pos;
    }

    PS::F64 getRSearch() const {
      return Parameter::search_rad;
    }
    
    void copyFromFP(const FPDPD& fp) {
      this->id          = fp.id;
      this->prop	= fp.prop;
      this->pos		= fp.pos;
    }
  };
};

namespace EPJ {
  struct DPD {
    PS::U32 id, prop;
    PS::F64vec pos, vel;
    std::array<PS::F64, Parameter::prop_num> dens;
    
    void copyFromFP(const FPDPD& fp) {
      this->id			= fp.id;
      this->prop		= fp.prop;
      this->pos			= fp.pos;
      this->vel			= fp.vel;
      this->dens = fp.density;
    }
    
    PS::F64vec getPos() const {
      return this->pos;
    }
    
    void setPos(const PS::F64vec& pos_) {
      this->pos = pos_;
    }
  };

  struct Density {
    PS::U32 id, prop;
    PS::F64vec pos;
    
    void copyFromFP(const FPDPD& fp) {
      this->id                  = fp.id;
      this->prop		= fp.prop;
      this->pos			= fp.pos;
    }
    
    PS::F64vec getPos() const {
      return this->pos;
    }
    
    void setPos(const PS::F64vec& pos_) {
      this->pos = pos_;
    }
  };
};
