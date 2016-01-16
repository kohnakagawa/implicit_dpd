#include <iostream>
#include <ctime>
#include <cstdlib>
#include "particle_simulator.hpp"

struct val {
  PS::U32 key;
  PS::U32 getKey() {
    return key;
  }
};

const int num = 100;
val v[num], v_buf[num];

int main() {
  srand((size_t)time(NULL));
  for(int i = 0; i < num; i++)
    v[i].key = rand();
  
  PS::RadixSort<PS::U32> rsorter;
  rsorter.lsdSort(v, v_buf, 0, num);

  for(int i = 0; i < num; i++)
    std::cout << v[i].key << std::endl;

  for(int i = 0; i < num - 1; i++)
    assert(v[i].key <= v[i + 1].key);
}
