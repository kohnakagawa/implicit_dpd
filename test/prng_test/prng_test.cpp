#include <iostream>
#include "particle_simulator.hpp"
#include "parameter.hpp"
#include "f_calculator.hpp"
#include <numeric>
#include <ctime>
#include <sys/time.h>
#include <iomanip>

constexpr PS::U32 sample_n = 100;
constexpr PS::U32 all_time = 100;

double time_diff(timeval& tv1, timeval& tv2) {
  return tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec) * 1.0e-6;
}

void gen_random_val(PS::F64* buffer, const PS::U32 m_seed) {
  PS::U32 cnt = 0;
  for(PS::U32 time = 0; time < all_time; time++) {
    for(PS::U32 i = 0; i < sample_n; i++) {
      for(PS::U32 j = i + 1; j < sample_n; j++) {
	Saru saru(i, j, time + m_seed);
	const PS::F64 buf = saru.nrml();
	buffer[cnt++] = buf;
      }
    }
  }
}

int main() {
  PS::U32 m_seed = static_cast<PS::U32>(time(NULL));

  PS::F64 *buffer = nullptr, *buffer_base = nullptr;
  const PS::U32 buf_size = sample_n * (sample_n - 1) / 2 * all_time;
  buffer = new PS::F64 [buf_size];
  buffer_base = new PS::F64 [buf_size];

  timeval tv1, tv2;
  gettimeofday(&tv1, NULL);
  gen_random_val(buffer, m_seed); //first
  gettimeofday(&tv2, NULL);

  std::cout << "generation time\n";
  std::cout << std::setprecision(15);
  std::cout << time_diff(tv1, tv2) << std::endl;;

  gen_random_val(buffer_base, m_seed); //genrate again
  
  for(PS::U32 i = 0; i < buf_size; i++)
    assert(buffer_base[i] == buffer[i]);
  
  //calculate mean and sd
  const PS::F64 mean = std::accumulate(buffer, buffer + buf_size, 0.0) / buf_size;
  PS::F64 sd = 0.0;
  for(PS::U32 i = 0; i < buf_size; i++)
    sd += (buffer[i] - mean) * (buffer[i] - mean);
  sd /= buf_size;
  sd = std::sqrt(sd);
  
  std::cout << "Mean / Stdevp\n";
  std::cout << mean << " " << sd << std::endl;
  
  std::ofstream fout("random_val.txt");
  for(PS::U32 i = 0; i < buf_size; i++)
    fout << buffer[i] << std::endl;

  delete [] buffer;
  delete [] buffer_base;
}
