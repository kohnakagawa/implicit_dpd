#pragma once

// #define CALC_HEIGHT // calculate membrane height or not
// #define PAIRWISE_DPD

#ifdef CHEM_MODE
#define LOCAL_CHEM_EVENT
#endif

#ifdef ENABLE_GPU_CUDA
#error "GPU mode is not available!"
enum {
  NUM_GPU_IN_ONE_NODE = 2,
};

#define USE_FLOAT_VEC
// #define USE_DOUBLE_VEC // Do not support now

#if defined (USE_DOUBLE_VEC) && defined (USE_FLOAT_VEC)
#error "Please choose USE_DOUBLE_VEC or USE_FLOAT_VEC."
#endif

#endif
