#include <iostream>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <sys/time.h>
#include <parallel/algorithm>
#include <cassert>
#include <tbb/parallel_sort.h>

double time_diff(timeval& tv1, timeval& tv2) {
  return tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec) * 1.0e-6;
}

typedef unsigned int uint;

template<typename Dtype>
void check(Dtype* array, const uint size) {
  for (uint i = 0; i < size - 1; i++)
    assert(array[i] <= array[i + 1]);
}

int main() {
  const int num = 1000000;
  uint *v0 = new uint [num];
  uint *v1 = new uint [num];
  
  srand((size_t)time(NULL));
  for(int i = 0; i < num; i++)
    v1[i] = v0[i] = (uint) rand();

  timeval tv1, tv2, tv3, tv4;
  
  gettimeofday(&tv1, nullptr);
  // __gnu_parallel::sort(v0, v0 + num);
  tbb::parallel_sort(v0, v0 + num);
  gettimeofday(&tv2, nullptr);

  printf("%f\n", time_diff(tv1, tv2));

  gettimeofday(&tv3, nullptr);
  std::sort(v1, v1 + num);
  gettimeofday(&tv4, nullptr);

  printf("%f\n", time_diff(tv3, tv4));
  
  check(v0, num);
  check(v1, num);

  delete [] v0;
  delete [] v1;
}
