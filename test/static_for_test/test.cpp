#include "../../static_for.hpp"
#include <iostream>
#include <cassert>
#include <ctime>

static constexpr int num = 10;

void clear(double* a,
	   double* b,
	   double* c)
{
  for(int i = 0; i < num; i++) {
    a[i] = i * 2;
    b[i] = i + 10;
    c[i] = i + 20;
  }
}

void check(const double* a,
	   const double* a_ref)
{
  for(int i = 0; i < num; i++)
    assert(a[i] == a_ref[i]);
}

double accumu1(const double* a,
	       const double* b)
{
  double sum_ref = 0.0;
  for(int i = 0; i < num; i++)
    sum_ref += a[i] * b[i];
  return sum_ref;
}

double accumu(const double* a,
	      const double* b)
{
  return static_for_accum<0, num, mul>(a, b);
}

int main() {
  double a[num], b[num], c[num];
  double a_ref[num], b_ref[num], c_ref[num];

  //addto
  clear(a, b, c);
  clear(a_ref, b_ref, c_ref);
  for(int i = 0; i < num; i++)
    a_ref[i] += b_ref[i];
  static_for<0, num, addto>(a, b);
  check(a, a_ref);
  
  //add
  clear(a, b, c);
  clear(a_ref, b_ref, c_ref);
  for(int i = 0; i < num; i++)
    c_ref[i] = a_ref[i] + b_ref[i];
  static_for<0, num, add>(c, a, b);
  check(c, c_ref);

  //accumulate
   clear(a, b, c);
  clear(a_ref, b_ref, c_ref);
  double sum, sum_ref = 0.0;
  
  clock_t start, end;
  
  start = clock();
  sum_ref = accumu1(a, b);
  end = clock();
  
  std::cout << end - start << std::endl;
  
  start = clock();  
  sum = accumu(a, b);
  end = clock();

  std::cout << end - start << std::endl;
  
  assert(sum == sum_ref);
  
  std::cout << "Success test\n.";
}
