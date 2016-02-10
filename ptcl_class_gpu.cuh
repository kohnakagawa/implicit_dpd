#pragma once

#include "parameter.hpp"

#ifdef USE_FLOAT_VEC

typedef float4 VecPos;
typedef float4 VecForce;
typedef float  Dtype;

#elif defined USE_DOUBLE_VEC
// Do not support now.
typedef double4 VecPos;
typedef double4 VecForce;
typedef double  Dtype;

#else

#error "Vector type is not specified in user_defs.h"

#endif

namespace RESULT {
  //NOTE: T should be double4 or float4.
  template <class T>
  struct ForceGPU {
    T acc;
    T press;
  };
  
  //NOTE: T should be  double or float.
  template <class T>
  struct DensityGPU {
    T dens[Parameter::prop_num];
  };
}

namespace EPI {
  template <class T>
  struct DPDGPU {
    T pos; PS::U32 id_;
    T vel; PS::U32 prop_;
    PS::S32 id_walk;
    double dens[Parameter::prop_num];
    
    __host__ __device__ __forceinline__ PS::U32& id() {
      return id_;
    }
    __host__ __device__ __forceinline__ const PS::U32& id() const {
      return id_;
    }

    __host__ __device__ __forceinline__ PS::U32& prop() {
      return prop_;
    }
    __host__ __device__ __forceinline__ const PS::U32& prop() const {
      return prop_;
    }
  };
  
  template <>
  struct DPDGPU<double4> {
    union {
      double4 pos;
      PS::U32 id_[8];
    };
    union {
      double4 vel;
      PS::U32 prop_[8];
    };
    //NOTE: id   == __double2hiuint(pos.w)
    //      prop == __double2hiuint(vel.w)
    
    PS::S32 id_walk;
    double dens[Parameter::prop_num];
    
    __host__ __device__ __forceinline__ PS::U32& id() {
      return id_[6];
    }
    __host__ __device__ __forceinline__ const PS::U32& id() const {
      return id_[6];
    }

    __host__ __device__ __forceinline__ PS::U32& prop() {
      return prop_[6];
    }
    __host__ __device__ __forceinline__ const PS::U32& prop() const {
      return prop_[6];
    }
  };

  template <>
  struct DPDGPU<float4> {
    union {
      float4 pos;
      PS::U32 id_[4];
    };
    union {
      float4 vel;
      PS::U32 prop_[4];
    };
    //NOTE: id   == __float_as_uint(pos.w)
    //      prop == __float_as_uint(vel.w)

    PS::S32 id_walk;
    float dens[Parameter::prop_num];

    __host__ __device__ __forceinline__ PS::U32& id() {
      return id_[3];
    }
    __host__ __device__ __forceinline__ const PS::U32& id() const {
      return id_[3];
    }

    __host__ __device__ __forceinline__ PS::U32& prop() {
      return prop_[3];
    }
    __host__ __device__ __forceinline__ const PS::U32& prop() const {
      return prop_[3];
    }
  };
  
  template <class T>
  struct DensityGPU {
    T pos; PS::U32 prop_;
    PS::S32 id_walk;

    __host__ __device__ __forceinline__ PS::U32& prop() {
      return prop_;
    }
    __host__ __device__ __forceinline__ const PS::U32& prop() const {
      return prop_;
    }
  };

  template <>
  struct DensityGPU<double4> {
    union {
      double4 pos;
      PS::U32 prop_[8];
    };
    //NOTE: prop == __double2hiuint(pos.w)
    PS::S32 id_walk;
    
    __host__ __device__ __forceinline__ PS::U32& prop() {
      return prop_[6];
    }
    __host__ __device__ __forceinline__ const PS::U32& prop() const {
      return prop_[6];
    }
  };
  
  template <>
  struct DensityGPU<float4> {
    union {
      float4 pos;
      PS::U32 prop_[4];
    };
    //NOTE: prop == __float_as_uint(pos.w)
    PS::S32 id_walk;
    
    __host__ __device__ __forceinline__ PS::U32& prop() {
      return prop_[3];
    }
    __host__ __device__ __forceinline__ const PS::U32& prop() const {
      return prop_[3];
    }
  };
} //End of namespace EPI

namespace EPJ {
  template <class T> //NOTE T should be double3
  struct DPDGPU {
    T pos; PS::U32 id_;
    T vel; PS::U32 prop_;
    double dens[Parameter::prop_num];
    
    __host__ __device__ __forceinline__ PS::U32& id() {
      return id_;
    }
    __host__ __device__ __forceinline__ const PS::U32& id() const {
      return id_;
    }

    __host__ __device__ __forceinline__ PS::U32& prop() {
      return prop_;
    }
    __host__ __device__ __forceinline__ const PS::U32& prop() const {
      return prop_;
    }
  };

  template <>
  struct DPDGPU<double4> {
    union {
      double4 pos;
      PS::U32 id_[8];
    };
    union {
      double4 vel;
      PS::U32 prop_[8];
    };
    //NOTE: id   == __double2hiuint(pos.w)
    //      prop == __double2hiuint(vel.w)
    double dens[Parameter::prop_num];
    
    __host__ __device__ __forceinline__ PS::U32& id() {
      return id_[6];
    }
    __host__ __device__ __forceinline__ const PS::U32& id() const {
      return id_[6];
    }

    __host__ __device__ __forceinline__ PS::U32& prop() {
      return prop_[6];
    }
    __host__ __device__ __forceinline__ const PS::U32& prop() const {
      return prop_[6];
    }
  };

  template <>
  struct DPDGPU<float4> {
    union {
      float4 pos;
      PS::U32 id_[4];
    };
    union {
      float4 vel;
      PS::U32 prop_[4];
    };    
    //NOTE: id   == __float_as_uint(pos.w)
    //      prop == __float_as_uint(vel.w)
    float dens[Parameter::prop_num];

    __host__ __device__ __forceinline__ PS::U32& id() {
      return id_[3];
    }
    __host__ __device__ __forceinline__ const PS::U32& id() const {
      return id_[3];
    }

    __host__ __device__ __forceinline__ PS::U32& prop() {
      return prop_[3];
    }
    __host__ __device__ __forceinline__ const PS::U32& prop() const {
      return prop_[3];
    }
  };
  
  template <class T>
  struct DensityGPU {
    T pos; PS::U32 prop_;

    __host__ __device__ __forceinline__ PS::U32& prop() {
      return prop_;
    }
    __host__ __device__ __forceinline__ const PS::U32& prop() const {
      return prop_;
    }
  };

  template <>
  struct DensityGPU<double4> {
    union {
      double4 pos;
      PS::U32 prop_[8];
    };
    //NOTE: prop == __double2hiuint(pos.w)
    
    __host__ __device__ __forceinline__ PS::U32& prop() {
      return prop_[6];
    }
    __host__ __device__ __forceinline__ const PS::U32& prop() const {
      return prop_[6];
    }
  };  
  
  template <>
  struct DensityGPU<float4> {
    union {
      float4 pos;
      PS::U32 prop_[4];
    };
    //NOTE: prop == __float_as_uint(pos.w)
    
    __host__ __device__ __forceinline__ PS::U32& prop() {
      return prop_[3];
    }
    __host__ __device__ __forceinline__ const PS::U32& prop() const {
      return prop_[3];
    }
  };  
} //end of namespace EPJ

