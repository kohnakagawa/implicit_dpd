#pragma once

template <int beg, int end, template<int i> class Func, typename... Args>
struct static_for_inner {
  static_assert(beg < end, "range error");
  void operator()(Args...args){
    Func<beg>()(args...);
    static_for_inner<beg+1, end, Func, Args...>()(args...);
  }
};

template <int end, template<int i> class Func, typename... Args>
struct static_for_inner<end, end, Func, Args...>{
  void operator()(Args...args){}
};

template <int beg, int end, template<int i> class Func, typename... Args>
void static_for(Args...args){
  static_for_inner<beg, end, Func, Args...>()(args...);
}

template<int i>
struct addto {
  template<typename T>
  void operator()(T* dst, const T* src) {
    dst[i] += src[i];
  }
};

template<int i>
struct add {
  template<typename T>
  void operator()(T* dst, const T* src0, const T* src1) {
    dst[i] = src0[i] + src1[i];
  }
};

template <int beg, int end, template<int i> class Func, typename... Args>
struct static_for_inner_accum {
  static_assert(beg < end, "range error");
  double operator()(Args...args){
    return Func<beg>()(args...) + static_for_inner_accum<beg+1, end, Func, Args...>()(args...);
  }
};

template <int end, template<int i> class Func, typename... Args>
struct static_for_inner_accum<end, end, Func, Args...>{
  double operator()(Args...args){ return 0.0; }
};

template <int beg, int end, template<int i> class Func, typename... Args>
double static_for_accum(Args...args){
  return static_for_inner_accum<beg, end, Func, Args...>()(args...);
}

template<int i>
struct mul {
  template<typename T>
  double operator()(const T* src0, const T* src1) {
    return src0[i] * src1[i];
  }
};
