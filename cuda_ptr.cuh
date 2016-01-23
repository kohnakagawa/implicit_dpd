#pragma once
 
#include <cstdio>
#include <cstdlib>
#include <thrust/fill.h>
#include <thrust/device_ptr.h>
#include "cuda_error.cuh"
 
template <typename T>
struct cuda_ptr{
  T *dev_ptr;
  T *host_ptr;
  int size;
  thrust::device_ptr<T> thrust_ptr;
  
  cuda_ptr(){
    dev_ptr  = nullptr;
    host_ptr = nullptr;
    size     = -1;
  }
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
