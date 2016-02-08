#include <iostream>
#include "particle_simulator.hpp"
#include "parameter.hpp"
#include "saruprngCUDA.cuh"
#include "cuda_ptr.cuh"
#include <numeric>
#include <ctime>
#include <sys/time.h>
#include <iomanip>

enum {
  THREADS = 70,
  BLOCKS  = 10,
};

constexpr PS::U32 sample_n = THREADS * BLOCKS;
constexpr PS::U32 sys_size = sample_n * sample_n;

double time_diff(timeval& tv1, timeval& tv2) {
  return tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec) * 1.0e-6;
}

__global__ void gen_random_val(float* buffer,
			       const PS::U32 m_seed)
{
  const uint tid = threadIdx.x + blockIdx.x * blockDim.x;
  for(uint i = 0; i < sample_n; i++) {
    //mi_id < mj_id
    const bool flag = (tid > i);
    const uint mi_id = flag ? i : tid;
    const uint mj_id = flag ? tid : i;
    SaruGPU saru(mi_id, mj_id, m_seed);
    const float buf = saru.nrml_f();
    
    const uint dst_id = tid * sample_n + i;
    buffer[dst_id] = buf;
  }
}

int main() {
  PS::U32 m_seed = static_cast<PS::U32>(time(NULL));
  
  cuda_ptr<float> buffer, buffer_base;
  buffer.allocate(sys_size);
  buffer_base.allocate(sys_size);

  buffer.set_val(10.0);
  buffer_base.set_val(10.0);
  
  timeval tv1, tv2;
  gettimeofday(&tv1, nullptr);
  gen_random_val<<<BLOCKS, THREADS>>>(buffer_base, m_seed);
  cudaDeviceSynchronize();
  gettimeofday(&tv2, nullptr);
  
  std::cout << "generation time\n";
  std::cout << std::setprecision(15);
  std::cout << time_diff(tv1, tv2) << std::endl;;
  
  gen_random_val<<< BLOCKS, THREADS>>>(buffer, m_seed);
  cudaDeviceSynchronize();  

  //copy from device to host
  buffer.dev2host(); buffer_base.dev2host();
  
  //check
  for(int i = 0; i < sys_size; i++)
    assert(buffer_base[i] == buffer[i]);
  
  std::ofstream fout("random_val.txt");
  fout << std::setprecision(8);
  for(PS::U32 i = 0; i < sample_n * sample_n; i++)
    fout << buffer[i] << std::endl;

  float mean = 0.0, sd = 0.0;
  for(int i = 0; i < sys_size; i++)
    mean += buffer[i];
  mean /= sys_size;

  for(int i = 0; i < sys_size; i++)
    sd += (buffer[i] - mean) * (buffer[i] - mean);
  sd /= sys_size;
  sd = std::sqrt(sd);
  
  std::cout << "samplen " << sys_size << std::endl;
  std::cout << mean << " " << sd << std::endl;

  buffer.deallocate();
  buffer_base.deallocate();
}