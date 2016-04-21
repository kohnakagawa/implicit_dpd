#pragma once
#include <iostream>

#define CHECK_EQ(val0, val1)						\
  do {									\
    if (val0 != val1) {							\
      std::cerr << "Is not equal.\n";					\
      std::cerr << "at " __FILE__ << " " << __LINE__ << std::endl;	\
      std::cerr << #val0 << " " << #val1 << std::endl;			\
      std::cerr << val0 << " " << val1 << std::endl;			\
      std::exit(1);							\
    }									\
  } while (false)

#define CHECK_LE(val0, val1)						\
  do {									\
    if (val0 > val1) {							\
      std::cerr << "val0 <= val1 is not statisfied.\n";			\
      std::cerr << "at " __FILE__ << " " << __LINE__ << std::endl;	\
      std::cerr << #val0 << " " << #val1 << std::endl;			\
      std::cerr << val0 << " " << val1 << std::endl;			\
      std::exit(1);							\
    }									\
  } while (false)

#define CHECK_LT(val0, val1)						\
  do {									\
    if (val0 >= val1) {							\
      std::cerr << "val0 < val1 is not statisfied.\n";			\
      std::cerr << "at " __FILE__ << " " << __LINE__ << std::endl;	\
      std::cerr << #val0 << " " << #val1 << std::endl;			\
      std::cerr << val0 << " " << val1 << std::endl;			\
      std::exit(1);							\
    }									\
  } while (false)

#define CHECK_GE(val0, val1)						\
  do {									\
    if (val0 < val1) {							\
      std::cerr << "val0 >= val1 is not statisfied.\n";			\
      std::cerr << "at " __FILE__ << " " << __LINE__ << std::endl;	\
      std::cerr << #val0 << " " << #val1 << std::endl;			\
      std::cerr << val0 << " " << val1 << std::endl;			\
      std::exit(1);							\
    }									\
  } while (false)

#define CHECK_GT(val0, val1)						\
  do {									\
    if (val0 <= val1) {							\
      std::cerr << "val0 > val1 is not statisfied.\n";			\
      std::cerr << "at " __FILE__ << " " << __LINE__ << std::endl;	\
      std::cerr << #val0 << " " << #val1 << std::endl;			\
      std::cerr << val0 << " " << val1 << std::endl;			\
      std::exit(1);							\
    }									\
  } while (false)

#define CHECK_FILE_OPEN(fin, name)				\
  do {								\
    if (!fin) {							\
      std::cerr << "Cannot open file " << name << std::endl;	\
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;	\
      std::exit(1);						\
    }								\
  } while (false)
