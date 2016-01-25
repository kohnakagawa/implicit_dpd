#include <stdio.h>

enum {
  Nx = 3,
  Ny = 4,
  Nz = 5,
};

texture<float, cudaTextureType2D, cudaReadModeElementType> texture2D;
texture<float, cudaTextureType3D, cudaReadModeElementType> texture3D;
texture<int2, cudaTextureType2D, cudaReadModeElementType> texture2D_d;

__global__ void kernel2D() {
  const int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid == 0) {
    for(int ix = 0; ix < Nx; ix++) {
      for(int iy = 0; iy < Ny; iy++) {
	const float val = tex2D(texture2D, iy, ix);
	const int2 v = tex2D(texture2D_d, iy, ix);
	const double val_d = __hiloint2double( v.y, v.x );
	printf("%f %f\n", val, val_d);
      }
    }
  }
}

__global__ void kernel3D() {
  const int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid == 0) {
    for(int ix = 0; ix < Nx; ix++) {
      for(int iy = 0; iy < Ny; iy++) {
	for(int iz = 0; iz < Nz; iz++) {
	  const float val = tex3D(texture3D, iz, iy, ix);
	  printf("%f\n", val);
	}
      }
    }
  }  
}

int main() {
  float host_2D[Nx][Ny]; //H W
  float host_3D[Nx][Ny][Nz]; //D H W
  double host_2D_d[Nx][Ny]; //H W

  int cnt = 0;
  for(int ix = 0; ix < Nx; ix++) {

  for(int iy = 0; iy < Ny; iy++) {
      host_2D[ix][iy] = cnt;
      host_2D_d[ix][iy] = cnt;
      cnt++;
    }
  }
  cnt = 0;
  for(int ix = 0; ix < Nx; ix++) {
    for(int iy = 0; iy < Ny; iy++) {
      for(int iz = 0; iz < Nz; iz++) {
        host_3D[ix][iy][iz] = cnt;
        cnt++;
      }
    }
  }
  
  cudaArray *cu_2D = nullptr, *cu_3D = nullptr, *cu_2D_d = nullptr;
  cudaChannelFormatDesc cdesc = cudaCreateChannelDesc<float>();
  cudaChannelFormatDesc cdesc_d = cudaCreateChannelDesc<int2>();
  
  cudaMallocArray(&cu_2D, &cdesc, Ny, Nx);
  cudaMalloc3DArray(&cu_3D, &cdesc, make_cudaExtent(Nz, Ny, Nx) );
  cudaMallocArray(&cu_2D_d, &cdesc_d, Ny, Nx);

  const size_t size2d = Nx * Ny * sizeof(float);
  cudaMemcpyToArray(cu_2D, 0, 0, host_2D, size2d, cudaMemcpyHostToDevice);

  const size_t size2d_d = Nx * Ny * sizeof(double);
  cudaMemcpyToArray(cu_2D_d, 0, 0, host_2D_d, size2d_d, cudaMemcpyHostToDevice);

  cudaMemcpy3DParms copyParams = {0};
  copyParams.srcPtr = make_cudaPitchedPtr(host_3D, Nz * sizeof(float), Nz, Ny);
  copyParams.dstArray = cu_3D;
  copyParams.extent = make_cudaExtent(Nz, Ny, Nx); //width height depth
  copyParams.kind   = cudaMemcpyHostToDevice;
  cudaMemcpy3D(&copyParams);

  texture2D.normalized = false;
  texture3D.normalized = false;
  
  cudaBindTextureToArray(texture2D, cu_2D, cdesc);
  cudaBindTextureToArray(texture3D, cu_3D, cdesc);
  cudaBindTextureToArray(texture2D_d, cu_2D_d, cdesc_d);
  
  kernel2D<<<32, 10>>>();
  kernel3D<<<32, 10>>>();

  cudaDeviceSynchronize();
  
  cudaFreeArray(cu_2D);
  cudaFreeArray(cu_3D);
}