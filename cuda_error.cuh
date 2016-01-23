#pragma once

template< typename T >
void check(T result, char const *const func, const char *const file, int const line){
  if(result){
    fprintf(stderr, "CUDA error at %s:%d code=%d \"%s\" \n",
	    file, line, static_cast<unsigned int>(result), func);
    cudaDeviceReset();
    // Make sure we call CUDA Device Reset before exiting
    exit(EXIT_FAILURE);
  }
}

#define checkCudaErrors(val) check((val), #val, __FILE__, __LINE__)

