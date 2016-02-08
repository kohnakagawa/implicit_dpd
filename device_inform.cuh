#pragma once

void print_device_inform(const int device_id) {
  cudaDeviceProp deviceProp;
  int driverVersion;
  int runtimeVersion;

  cudaGetDeviceProperties(&deviceProp, device_id);
  cudaDriverGetVersion(&driverVersion);
  cudaRuntimeGetVersion(&runtimeVersion);

  const float peak_mem_bw = (float) deviceProp.memoryClockRate * 1e3 * 2 * deviceProp.memoryBusWidth / 8;

  printf("\n");
  printf("******** GPU Information ********\n");
  printf("CUDA Driver Version, %d.%d\n", driverVersion/1000, (driverVersion%100)/10 );
  printf("CUDA Runtime Version, %d.%d\n", runtimeVersion/1000, (runtimeVersion%100)/10);
  printf("CUDA Capability Major/Minor version number, %d.%d\n", deviceProp.major, deviceProp.minor);
  printf( "GPU - name, %s\n", deviceProp.name );
  printf( "GPU - # SMs, %d\n", deviceProp.multiProcessorCount );
  printf( "GPU - core freq (GHz), %lf\n", deviceProp.clockRate * 1e-6 );
  printf( "GPU - L2 size (MiB), %lf\n", (float) deviceProp.l2CacheSize / 1024 / 1024 );
  printf( "GPU - memory size (GiB), %lf\n", (float) deviceProp.totalGlobalMem / 1024 / 1024 / 1024 );
  printf( "GPU - memory bus width (bit), %d\n", deviceProp.memoryBusWidth );
  printf( "GPU - memory bus freq (GHz), %lf\n", (float) deviceProp.memoryClockRate * 1e-6 );
  printf( "GPU - theoretical peak memory bandwidth (GB/s), %lf\n", peak_mem_bw * 1e-9 );
  printf( "GPU - ECC, %s\n", (deviceProp.ECCEnabled)?"On":"Off" );
  printf( "\n" );
}