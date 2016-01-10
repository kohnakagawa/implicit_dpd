#pragma once

#include <array>

namespace RESULT {
  struct ForceDPDDensity {
    PS::F64vec acc, press;
    std::array<PS::F64, Parameter::prop_num> dens;
    
    void clear(){
      acc   = 0.0;
      press = 0.0;
      dens.fill(0.0);
    }
  };
};

#include "full_ptcl.hpp"

namespace EPI {
  struct DPD_Density {
    PS::U32 id, id_loc, prop;
    PS::F64vec pos, vel;
    std::array<PS::F64, Parameter::prop_num> dens;

    PS::F64vec getPos() const {
      return this->pos;
    }

    PS::F64 getRSearch() const {
      return Parameter::search_rad;
    }

    void copyFromFP(const FPDPD& fp) {
      this->id		= fp.id;
      this->id_loc      = fp.id_loc;
      this->prop	= fp.prop;
      this->pos		= fp.pos;
      this->vel		= fp.vel;
      this->dens        = fp.density;
    }

  };
};

namespace EPJ {
  struct DPD_Density {
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
};
