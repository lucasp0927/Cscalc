#include "cumtrapz.h"

int inline index(int i,int j,int k,int N)
{
  //  return 2*k+N*2*num*j+2*num*i;
  return k*N*N*2+N*2*j+2*i;
}
void addorder(double* T,double* D,double* env, double* A,int N, int num,double dt)
{
    typedef struct{ double re; double im; } complex16;
    double tmp[N*N*2];
    double result[N*N*2];
    complex16 alpha;
    alpha.re = 1.0;
    alpha.im = 0.0;
    complex16 beta;
    beta.re = 0.0;
    beta.im = 0.0;
    complex16 ef;
    for (int i = 0; i < num; ++i)
    {
	cblas_zcopy(N*N,T,1,tmp,1);
	ef.re = env[i];
	ef.im = 0.0;
	cblas_zaxpy(N*N,&ef,D,1,tmp,1);
	cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,N,N,N,&alpha,tmp,N,&A[index(0,0,i,N)],N,&beta,result,N);
	cblas_zcopy(N*N,result,1,&A[index(0,0,i,N)],1);
    }
    cumtrapz(A,dt,N,num);
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
