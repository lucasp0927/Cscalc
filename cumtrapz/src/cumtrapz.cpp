#include "cumtrapz.h"


void cumtrapz(double* A,double dt,int N,int num)
{
// # pragma omp parallel for
//   for (int i = 0; i < 2*N*N; ++i)
//     {
//       result[i] = A[i]+B[i];
//     }
  // real part
  double tmp[num];
  for (int i = 0; i < N; ++i)
    {
      for (int j = 0; j < N; ++j)
        {
          if (i==j)
            tmp[0] = 1.0;
          else
            tmp[0] = 0.0;
          for (int k = 1; k < num; ++k)
            tmp[k] = tmp[k-1] + (A[(k*N*N+N*j+i)*2]+A[(N*N*(k-1)+N*j+i)*2])*dt/2.0;
          for (int k = 0; k < num; ++k)
            A[(k*N*N+N*j+i)*2] =tmp[k];     
        }
    }
  // imag part
  for (int i = 0; i < N; ++i)
    {
      for (int j = 0; j < N; ++j)
        {
          tmp[0] = 0.0;
          for (int k = 1; k < num; ++k)
            tmp[k] = tmp[k-1] + (A[(k*N*N+N*j+i)*2+1]+A[((k-1)*N*N+N*j+i)*2+1])*dt/2.0;
          for (int k = 0; k < num; ++k)
            A[(k*N*N+N*j+i)*2+1] =tmp[k];     
        }
    }  
}




