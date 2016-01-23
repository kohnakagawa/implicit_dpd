#pragma once

#define CALC_HEIGHT // calculate membrane height or not
#define CHEM_MODE   // chem mode on or off
#define USE_GPU     // use gpu or not

#ifdef USE_GPU

#define USE_FLOAT_VEC
//#define USE_DOUBLE_VEC

#if defined (USE_DOUBLE_VEC) && defined (USE_FLOAT_VEC)
#error "Please choose USE_DOUBLE_VEC or USE_FLOAT_VEC."
#endif

#endif
