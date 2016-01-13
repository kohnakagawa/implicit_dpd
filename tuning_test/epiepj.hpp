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

#include "full_ptcl.hpp"

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
