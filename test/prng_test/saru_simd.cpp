#include <iostream>
#include <array>
#include <vector>
#include "saruprng.hpp"
#include "saruprngAVX2.hpp"

std::array<int, 3> get_position(int key, const int dim) {
  const auto k = key % dim;
  auto tmp = key / dim;
  const auto j = tmp % dim;
  const auto i = tmp / dim;
  return {i, j, k};
}

void print256d(v4df v) {
  union {
    v4df v;
    double elem[4];
  } tmp;
  tmp.v = v;
  for (int i = 0; i < 4; i++) {
    std::cout << tmp.elem[i] << std::endl;
  }
}

void function() {
  v4df v0 = _mm256_set_pd(0.4, 0.3, 0.2, 0.1);
  v4df v1 = _mm256_set_pd(0.2, 0.3, 0.4, 0.5);

  v4df in_cos  = (v4df)_mm256_set1_pd(2.0 * M_PI) * v0;
  v4df in_sqrt = (v4df)_mm256_set1_pd(-2.0) * (v4df)_mm256_log_pd(v1);
  v4df ret = (v4df)_mm256_sqrt_pd(in_sqrt) * (v4df)_mm256_cos_pd(in_cos);
  print256d(ret);

  std::cout << std::endl;

  const double in_cs = 2.0 * M_PI * 0.2;
  const double in_sq = -2.0 * std::log(0.4);
  std::cout << std::sqrt(in_sq) * std::cos(in_cs) << std::endl;
}

int main(int argc, char* argv[]) {
  int num_test = 100;
  if (argc == 2) {
    num_test = std::atoi(argv[1]);
  }

  const auto array_size = num_test * num_test * num_test;
  std::cerr << "# of test cases is " << array_size << std::endl;
  std::vector<double> array(array_size), array_ref(array_size);

#if 0
  for (int i = 2; i < 6; i++) {
    Saru saru(0, 1, i);
    std::cout << saru.nrml() << std::endl;
  }

  std::cout << std::endl;

  v8si vs1 = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 0, 0);
  v8si vs2 = _mm256_set_epi32(0, 0, 0, 0, 1, 1, 1, 1);
  v8si vs3 = _mm256_set_epi32(0, 0, 0, 0, 5, 4, 3, 2);
  SaruAVX2 saru(vs1, vs2, vs3);
  print256d(saru.vnrml());
#else
  // reference
  int cnt = 0;
  for (int i = 0; i < num_test; i++) {
    for (int j = 0; j < num_test; j++) {
      for (int k = 0; k < num_test; k++) {
        Saru saru(i, j, k);
        array_ref[cnt++] = saru.nrml();
      }
    }
  }

  // simd
  const int num_loop = num_test * num_test * num_test;
  for (cnt = 0; cnt < num_loop; cnt += 4) {
    const auto pos0 = get_position(cnt    , num_test);
    const auto pos1 = get_position(cnt + 1, num_test);
    const auto pos2 = get_position(cnt + 2, num_test);
    const auto pos3 = get_position(cnt + 3, num_test);

    v8si vs1 = _mm256_set_epi32(0, 0, 0, 0, pos3[0], pos2[0], pos1[0], pos0[0]);
    v8si vs2 = _mm256_set_epi32(0, 0, 0, 0, pos3[1], pos2[1], pos1[1], pos0[1]);
    v8si vs3 = _mm256_set_epi32(0, 0, 0, 0, pos3[2], pos2[2], pos1[2], pos0[2]);
    SaruAVX2 saru(vs1, vs2, vs3);

    v4df vec = saru.vnrml();
    _mm256_storeu_pd(&array[cnt], vec);
  }

  std::cout << "(array - array_ref)/array_ref \n";
  for (int i = 0; i < 30; i++) {
    std::cout << (array[i] - array_ref[i]) / array_ref[i] << std::endl;
  }
#endif
}
