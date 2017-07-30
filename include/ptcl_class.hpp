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
    PS::F64vec nei_pos_sum[2];
    PS::U32 nei_cnt[2];

    void clear() {
      dens.fill(0.0);
      nei_pos_sum[0] = nei_pos_sum[1] = 0.0;
      nei_cnt[0] = nei_cnt[1] = 0;
    }
  };
};

struct FPDPD {
  PS::U32 id, prop, amp_id, unit;
  PS::F64vec pos, delta_sumr;
  PS::F64vec vel, vel_buf;
  PS::F64vec acc;
  PS::F64vec press;
  PS::F64vec nei_cm_pos[2];
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
    nei_cm_pos[Parameter::Hyphil] = dens.nei_pos_sum[Parameter::Hyphil] / dens.nei_cnt[Parameter::Hyphil];
    nei_cm_pos[Parameter::Hyphob] = dens.nei_pos_sum[Parameter::Hyphob] / dens.nei_cnt[Parameter::Hyphob];
  }

  //for I/O
  void readAscii(FILE *fp) {
    char buf;
    assert(0 != fscanf(fp, "%c %lf %lf %lf %u %u %u %u %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
		       &buf, &(pos.x), &(pos.y), &(pos.z),
		       &id, &prop, &amp_id, &unit,
		       &(vel.x), &(vel.y), &(vel.z), &(vel_buf.x), &(vel_buf.y), &(vel_buf.z),
		       &(acc.x), &(acc.y), &(acc.z)));
    delta_sumr = press = nei_cm_pos[0] = nei_cm_pos[1] = 0.0;
    density.fill(0.0);
  }
  void writeAscii(FILE *fp) const {
    fprintf(fp, "%c %.15g %.15g %.15g %u %u %u %u %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n",
	    Parameter::atom_type[prop], pos.x, pos.y, pos.z,
	    id, prop, amp_id, unit,
	    vel.x, vel.y, vel.z, vel_buf.x, vel_buf.y, vel_buf.z,
	    acc.x, acc.y, acc.z);
  }
  void readFromString(const char* line) {
    char buf;
    assert(0 != sscanf(line, "%c %lf %lf %lf %u %u %u %u %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
		       &buf, &(pos.x), &(pos.y), &(pos.z),
		       &id, &prop, &amp_id, &unit,
		       &(vel.x), &(vel.y), &(vel.z), &(vel_buf.x), &(vel_buf.y), &(vel_buf.z),
		       &(acc.x), &(acc.y), &(acc.z)));
    delta_sumr = press = nei_cm_pos[0] = nei_cm_pos[1] = 0.0;
    density.fill(0.0);
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
    PS::U32 prop;
    PS::F64vec pos;

    PS::F64vec getPos() const {
      return this->pos;
    }

    PS::F64 getRSearch() const {
      return Parameter::search_rad;
    }

    void copyFromFP(const FPDPD& fp) {
      this->prop	= fp.prop;
      this->pos		= fp.pos;
    }
  };
};

namespace EPJ {
  struct DPD {
    PS::U32 id, prop, amp_id, unit;
    PS::F64vec pos, vel;
    std::array<PS::F64, Parameter::prop_num> dens;

    void copyFromFP(const FPDPD& fp) {
      this->id			= fp.id;
      this->prop		= fp.prop;
      this->amp_id		= fp.amp_id;
      this->unit		= fp.unit;
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
    PS::U32 prop;
    PS::F64vec pos;

    void copyFromFP(const FPDPD& fp) {
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
