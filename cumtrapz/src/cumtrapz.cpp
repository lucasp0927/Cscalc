#include "cumtrapz.h"

int inline index(int i,int j,int k,int N,int num)
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
          tmp[1] = tmp[0] + (A[index(i,j,1,N,num)]+A[index(i,j,0,N,num)])*dt/2.0;          
          for (int k = 2; k < num; ++k)
            tmp[k] = tmp[k-2] + (A[index(i,j,k,N,num)]+A[index(i,j,k-1,N,num)]*4+A[index(i,j,k-2,N,num)])*dt/3.0;
          for (int k = 0; k < num; ++k)
            A[index(i,j,k,N,num)] =tmp[k];
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
          tmp[1] = tmp[0] + (A[index(i,j,1,N,num)+1]+A[index(i,j,0,N,num)+1])*dt/2.0;          
          for (int k = 2; k < num; ++k)
            tmp[k] = tmp[k-2] + (A[index(i,j,k,N,num)+1]+A[index(i,j,k-1,N,num)+1]*4+A[index(i,j,k-2,N,num)+1])*dt/3.0;          
          for (int k = 0; k < num; ++k)
            A[index(i,j,k,N,num)+1] =tmp[k];
        }
    }
}
