#pragma once

//#define CALC_HEIGHT // calculate membrane height or not
//#define CHEM_MODE   // chem mode on or off
//#define PAIRWISE_DPD

#ifdef ENABLE_GPU_CUDA

enum {
  NUM_GPU_IN_ONE_NODE = 2,
};

#define USE_TEXTURE_MEM

#define USE_FLOAT_VEC
// #define USE_DOUBLE_VEC // Do not support now

#if defined (USE_DOUBLE_VEC) && defined (USE_FLOAT_VEC)
#error "Please choose USE_DOUBLE_VEC or USE_FLOAT_VEC."
#endif

#endif
