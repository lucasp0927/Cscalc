#include "matrixMul.h"
// void deviceVerify ()
// {
//   int devID;
//   cudaDeviceProp props;
//   // get number of SMs on this GPU
//   cudaGetDevice(&devID);
//   cudaGetDeviceProperties(&props, devID);
//   printf("Device %d: \"%s\" with Compute %d.%d capability\n", devID, props.name, props.major, props.minor);
// }

void zdot(double* A,double* B,int N,double* result)
{
  unsigned int mem_size = sizeof(cuDoubleComplex)*N*N;
  cuDoubleComplex* d_result,* d_A,* d_B;
  cudaMalloc((void**) &d_result,mem_size);
  cudaMalloc((void**) &d_A,mem_size);
  cudaMalloc((void**) &d_B,mem_size);
  cudaMemcpy(d_A, A, mem_size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_B, B, mem_size, cudaMemcpyHostToDevice);  
  const cuDoubleComplex alpha = make_cuDoubleComplex(1.0,0.0);
  const cuDoubleComplex beta = make_cuDoubleComplex(0.0,0.0);
  cudaEventRecord(start, 0);

  cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N, &alpha, d_B, N, d_A, N, &beta, d_result, N);
  cudaMemcpy(result, d_result, mem_size, cudaMemcpyDeviceToHost);
  cudaFree(d_A);
  cudaFree(d_B);
  cudaFree(d_result);  
}

// void cdot(float* A,float* B,int N,float* result)
// {
//   cudaEvent_t start, stop;
//   float time;
//   cudaEventCreate(&start);
//   cudaEventCreate(&stop);
  
//   cublasHandle_t handle;
//   cublasCreate(&handle);
// //checkError(cublasCreate(&handle), "cublasCreate() error!\n");  
//   unsigned int mem_size = sizeof(cuComplex)*N*N;
  
//   cuComplex* d_result,* d_A,* d_B;
//   cudaMalloc((void**) &d_result,mem_size);
//   cudaMalloc((void**) &d_A,mem_size);
//   cudaMalloc((void**) &d_B,mem_size);
//   cudaMemcpy(d_A, A, mem_size, cudaMemcpyHostToDevice);
//   cudaMemcpy(d_B, B, mem_size, cudaMemcpyHostToDevice);  
//   const cuComplex alpha = make_cuComplex(1.0,0.0);
//   const cuComplex beta = make_cuComplex(0.0,0.0);
//   cudaEventRecord(start, 0);

//   cublasCgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N, &alpha, d_B, N, d_A, N, &beta, d_result, N);
//   cudaEventRecord(stop, 0);
//   cudaEventSynchronize(stop);
//   cudaEventElapsedTime(&time, start, stop);
//   printf ("Time for the cublasZgemm: %f ms\n", time);
  
//   cudaMemcpy(result, d_result, mem_size, cudaMemcpyDeviceToHost);
  
//   cudaFree(d_A);
//   cudaFree(d_B);
//   cudaFree(d_result);  
// }

