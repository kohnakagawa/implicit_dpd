#pragma once

#include "parameter.hpp"

//NOTE: T should be float3 or double3 or float4
namespace RESULT {
  template <class T>
  struct ForceGPU {
    T acc;
    T press;
  };
  
  template <class T>
  struct DensityGPU {
    T dens[Parameter::prop_num];
  };
}

namespace EPI {
  template <class T> //NOTE T should be double3
  struct DPDGPU {
    PS::U32 id_, prop_;
    PS::S32 id_walk;
    T pos, vel;
    T dens[Parameter::prop_num];
    
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
  struct DPDGPU<float4> {
    union {
      float4 pos;
      PS::U32 id_[4];
    };
    union {
      float4 vel;
      PS::U32 prop_[4];
    };
    //NOTE: id   == __float_as_int(pos.w)
    //      prop == __float_as_int(vel.w)

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
    PS::U32 prop_;
    PS::S32 id_walk;
    T pos;

    __host__ __device__ __forceinline__ PS::U32& prop() {
      return prop_;
    }
    __host__ __device__ __forceinline__ const PS::U32& prop() const {
      return prop_;
    }
  };
  template <>
  struct DensityGPU<float4> {
    union {
      float4 pos;
      PS::U32 prop_[4];
    };
    //NOTE: prop == __float_as_int(pos.w)
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
  template <class T> //NOTE T should be float3 or double3
  struct DPDGPU {
    PS::U32 id_, prop_;
    T pos, vel;
    T dens[Parameter::prop_num];
    
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
  struct DPDGPU<float4> {
    union {
      float4 pos;
      PS::U32 id_[4];
    };
    union {
      float4 vel;
      PS::U32 prop_[4];
    };    
    //NOTE: id   == __float_as_int(pos.w)
    //      prop == __float_as_int(vel.w)
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
    PS::U32 prop_;
    T pos;

    __host__ __device__ __forceinline__ PS::U32& prop() {
      return prop_;
    }
    __host__ __device__ __forceinline__ const PS::U32& prop() const {
      return prop_;
    }
  };
  template <>
  struct DensityGPU<float4> {
    union {
      float4 pos;
      PS::U32 prop_[4];
    };
    //NOTE: prop == __float_as_int(pos.w)
    
    __host__ __device__ __forceinline__ PS::U32& prop() {
      return prop_[3];
    }
    __host__ __device__ __forceinline__ const PS::U32& prop() const {
      return prop_[3];
    }
  };  
} //end of namespace EPJ

#ifdef USE_FLOAT_VEC

typedef float4 VecPos;
typedef float4 VecForce;
typedef float  Dtype;

#elif defined USE_DOUBLE_VEC

typedef double3 VecPos;
typedef double3 VecForce;
typedef double  Dtype;

#else

#error "Vector type is not specified in user_defs.h"

#endif
