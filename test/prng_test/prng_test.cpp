#include <iostream>
#include "particle_simulator.hpp"
#include "parameter.hpp"
#include "f_calculator.hpp"
#include <numeric>
#include <ctime>

constexpr PS::U32 sample_n = 200;
constexpr PS::U32 all_time = 200;

void gen_random_val(PS::F64* buffer, const PS::U32 m_seed) {
  PS::U32 cnt = 0;
  for(PS::U32 time = 0; time < all_time; time++) {
    for(PS::U32 i = 0; i < sample_n; i++) {
      for(PS::U32 j = i + 1; j < sample_n; j++) {
	Saru saru(i, j, time + m_seed);
	const PS::F64 buf = std::sqrt(-2.0f * std::log(saru.f() ) ) * std::cos(2.0f * M_PI * saru.f());
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

  gen_random_val(buffer, m_seed); //first
  gen_random_val(buffer_base, m_seed); //again genrate
  
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
