#pragma once

namespace alloc_util {
  //NOTE: T should not be pointer type.
  // ptr[Nx][Ny]
  template<class T>
  inline void allocate2D(const int Nx, const int Ny, T** ptr) {
    if(ptr == nullptr) {
      allocate2D_nocheck(Nx, Ny, ptr);
    } else {
      deallocate2D(ptr);
      allocate2D_nocheck(Nx, Ny, ptr);
    }
  }

  template<class T>
  inline void allocate2D_nocheck(const int Nx, const int Ny, T** ptr) {
    ptr = new T* [Nx];
    ptr[0] = new T[Nx * Ny];
    for(int i = 1; i < Nx; i++)
      ptr[i] = ptr[0] + i * Ny;
  }
  
  template<class T>
  inline void deallocate2D(T** ptr) {
    delete [] ptr[0];
    delete [] ptr;
    ptr = nullptr;
  }
};
