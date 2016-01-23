#pragma once

//NOTE: T should be float3 or double3
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
};

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
    float4 pos, vel; //NOTE: id   == __float_as_int(pos.w)
                     //      prop == __float_as_int(vel.w)
    PS::S32 id_walk;
    float dens[Parameter::prop_num];

    __host__ __device__ __forceinline__ PS::U32& id() {
      return *reinterpret_cast<PS::U32>(&(pos.w));
    }
    __host__ __device__ __forceinline__ const PS::U32& id() const {
      return *reinterpret_cast<PS::U32>(&(pos.w));
    }

    __host__ __device__ __forceinline__ PS::U32& prop() {
      return *reinterpret_cast<PS::U32>(&(vel.w));
    }
    __host__ __device__ __forceinline__ const PS::U32& prop() const {
      return *reinterpret_cast<PS::U32>(&(vel.w));
    }
  }
  
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
    float4 pos; //NOTE: prop == __float_as_int(pos.w)
    
    __host__ __device__ __forceinline__ PS::U32& prop() {
      return *reinterpret_cast<PS::U32>(&(pos.w));
    }
    __host__ __device__ __forceinline__ const PS::U32& prop() const {
      return *reinterpret_cast<PS::U32>(&(pos.w));
    }
  };
};

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
    float4 pos, vel; //NOTE: id   == __float_as_int(pos.w)
                     //      prop == __float_as_int(vel.w)
    float dens[Parameter::prop_num];

    __host__ __device__ __forceinline__ PS::U32& id() {
      return *reinterpret_cast<PS::U32>(&(pos.w));
    }
    __host__ __device__ __forceinline__ const PS::U32& id() const {
      return *reinterpret_cast<PS::U32>(&(pos.w));
    }

    __host__ __device__ __forceinline__ PS::U32& prop() {
      return *reinterpret_cast<PS::U32>(&(vel.w));
    }
    __host__ __device__ __forceinline__ const PS::U32& prop() const {
      return *reinterpret_cast<PS::U32>(&(vel.w));
    }
  }
  
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
    float4 pos; //NOTE: prop == __float_as_int(pos.w)
    
    __host__ __device__ __forceinline__ PS::U32& prop() {
      return *reinterpret_cast<PS::U32>(&(pos.w));
    }
    __host__ __device__ __forceinline__ const PS::U32& prop() const {
      return *reinterpret_cast<PS::U32>(&(pos.w));
    }
  };  
};

#ifdef USE_FLOAT_VEC

typedef float4 VecPos;
typedef float3 VecForce;

#elif defined USE_DOUBLE_VEC

typedef double3 VecPos;
typedef double3 VecForce;

#else

#error "Vector type is not specified in user_defs.h"

#endif


