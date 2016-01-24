#pragma once
 
#include <cstdio>
#include <cstdlib>
#include <thrust/fill.h>
#include <thrust/device_ptr.h>
#include "cuda_error.cuh"
 
template <typename T>
struct cuda_ptr{
  T *dev_ptr = nullptr;
  T *host_ptr = nullptr;
  int size = -1;
  thrust::device_ptr<T> thrust_ptr;
  
  cuda_ptr(){}
  ~cuda_ptr(){ deallocate(); }
 
  void allocate(const int size_){
    size = size_;
    checkCudaErrors(cudaMalloc((void**)&dev_ptr, size * sizeof(T)));
    checkCudaErrors(cudaMallocHost((void**)&host_ptr, size * sizeof(T)));
    thrust_ptr = thrust::device_pointer_cast(dev_ptr);
  }
  
  void host2dev(const int beg, const int count) {
    checkCudaErrors(cudaMemcpy(dev_ptr + beg, host_ptr + beg, count * sizeof(T), cudaMemcpyHostToDevice));
  }
  void host2dev() {this->host2dev(0, size);}
  void host2dev_async(const int beg, const int count, cudaStream_t& strm){
    checkCudaErrors(cudaMemcpyAsync(dev_ptr + beg, host_ptr + beg, count * sizeof(T), cudaMemcpyHostToDevice, strm));
  }
  
  void dev2host(const int beg, const int count){
    checkCudaErrors(cudaMemcpy(host_ptr + beg, dev_ptr + beg, count * sizeof(T), cudaMemcpyDeviceToHost));
  }
  void dev2host() {this->dev2host(0, size);}
  void dev2host_async(const int beg, const int count, cudaStream_t& strm){
    checkCudaErrors(cudaMemcpyAsync(host_ptr + beg, dev_ptr + beg, count * sizeof(T), cudaMemcpyDeviceToHost, strm));
  }
 
  void set_val(const T val){
    for(int i = 0; i < size; i++) host_ptr[i] = val;
    thrust::fill(thrust_ptr, thrust_ptr + size, val);
  }
  void set_val(const int beg, const int count, const T val){
    const int end = beg + count;
    for(int i = beg; i < end; i++) host_ptr[i] = val;
    thrust::device_ptr<T> beg_ptr = thrust_ptr + beg;
    thrust::fill(beg_ptr, beg_ptr + count, val);
  }
  
  const T & operator [] (const int i) const {
    return host_ptr[i];
  }
  
  T &operator [] (const int i) {
    return host_ptr[i];
  }
  
  operator T* (){
    return dev_ptr;
  }

private:
  void deallocate(){
    checkCudaErrors(cudaFree(dev_ptr));
    checkCudaErrors(cudaFreeHost(host_ptr));
  }
};

//for texture memory
template<typename T, size_t Nx, size_t Ny> //Height Width
struct cuda2D_ptr {
  size_t size[2] = {Nx, Ny};
  T host[Nx][Ny];
  cudaArray *dev_ptr = nullptr;
  cudaChannelFormatDesc cdesc;
  
  cuda2D_ptr(cudaChannelFormatDesc cdesc_) {
    cdesc = cdesc_;
    allocate();
  }
  ~cuda2D_ptr() {
    deallocate();
  }
  
  template<typename U>
  void host2host(U cpy[Nx][Ny]) {
    for(size_t i = 0; i < Nx; i++)
      for(size_t j = 0; j < Ny; j++)
	host[i][j] = (T)cpy[i][j];
  }
  
  void host2dev() {
    cudaMemcpyToArray(dev_ptr, 0, 0, host, Nx * Ny * sizeof(T), cudaMemcpyHostToDevice);
  }
  
private:
  void allocate() {
    cudaMallocArray(&dev_ptr, &cdesc, Ny, Nx);
  }

  void deallocate() {
    cudaFreeArray(dev_ptr);
  }
};

template<typename T, size_t Nx, size_t Ny, size_t Nz>
struct cuda3D_ptr {
  size_t size[3] = {Nx, Ny, Nz};
  T host[Nx][Ny][Nz];
  cudaArray *dev_ptr = nullptr;
  cudaChannelFormatDesc cdesc;
  cudaExtent ext = {Nz, Ny, Nx};
  
  cuda3D_ptr(cudaChannelFormatDesc cdesc_) {
    cdesc = cdesc_;
    allocate();
  }
  ~cuda3D_ptr() {
    deallocate();
  }

  template<typename U>
  void host2host(U cpy[Nx][Ny][Nz]) {
    for(size_t i = 0; i < Nx; i++)
      for(size_t j = 0; j < Ny; j++)
	for(size_t k = 0; k < Nz; k++)
	  host[i][j][k] = (T)cpy[i][j][k];
  }

  void host2dev() {
    cudaMemcpy3DParms copyParams = {0};
    copyParams.srcPtr = make_cudaPitchedPtr(host, Nz * sizeof(T), Nz, Ny);
    copyParams.dstArray = dev_ptr;
    copyParams.extent = ext;
    copyParams.kind   = cudaMemcpyHostToDevice;
    cudaMemcpy3D(&copyParams);
  }

private:
  void allocate() {
    cudaMalloc3DArray(&dev_ptr, &cdesc, ext);
  }

  void deallocate() {
    cudaFreeArray(dev_ptr);
  }
};