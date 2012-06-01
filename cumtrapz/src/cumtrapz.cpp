#include "cumtrapz.h"

int inline index(int i,int j,int k,int N)
{
  //  return 2*k+N*2*num*j+2*num*i;
  return k*N*N*2+N*2*j+2*i;
}

void cumtrapz(double* A,double dt,int N,int num)
{
  // real part
  double tmp[num];
  int i;
  int j;
# pragma omp parallel for  \
  private(i,j,tmp)
  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          if (i==j)
            tmp[0] = 1.0;
          else
            tmp[0] = 0.0;
          tmp[1] = tmp[0] + (A[index(i,j,1,N)] + A[index(i,j,0,N)])*dt/2.0;          
          for (int k = 2; k < num; ++k)
            tmp[k] = tmp[k-2] + (A[index(i,j,k,N)] + A[index(i,j,k-1,N)]*4 + A[index(i,j,k-2,N)])*dt/3.0;
          for (int k = 0; k < num; ++k)
            A[index(i,j,k,N)] =tmp[k];
        }
    }
  // imag part
# pragma omp parallel for  \
  private(i,j,tmp)
  for (i = 0; i < N; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          tmp[0] = 0.0;
          tmp[1] = tmp[0] + (A[index(i,j,1,N)+1]+A[index(i,j,0,N)+1])*dt/2.0;          
          for (int k = 2; k < num; ++k)
            tmp[k] = tmp[k-2] + (A[index(i,j,k,N)+1]+A[index(i,j,k-1,N)+1]*4+A[index(i,j,k-2,N)+1])*dt/3.0;          
          for (int k = 0; k < num; ++k)
            A[index(i,j,k,N)+1] =tmp[k];
        }
    }
}
