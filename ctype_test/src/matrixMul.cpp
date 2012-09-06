#include "matrixMul.h"


void zdot(double* A,double* B,int N,double* result)
{
  // unsigned int mem_size = sizeof(cuDoubleComplex)*N*N;
  // cuDoubleComplex* d_result,* d_A,* d_B;
# pragma omp parallel for
  for (int i = 0; i < 2*N*N; ++i)
    {
      result[i] = A[i]+B[i];
    }
  // cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N, &alpha, d_B, N, d_A, N, &beta, d_result, N);
  // cudaMemcpy(result, d_result, mem_size, cudaMemcpyDeviceToHost);
}




