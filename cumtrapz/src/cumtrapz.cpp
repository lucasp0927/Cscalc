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
# pragma omp parallel for  \
  private(tmp)
   for (int i = 0; i < N; ++i)
    {
      for (int j = 0; j < N; ++j)
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
  private(tmp)
  for (int i = 0; i < N; ++i)
    {
      for (int j = 0; j < N; ++j)
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
/*
char *int2bin(int a, char *buffer, int buf_size, int& size) {
    buffer += (buf_size - 1);

    for (int i = (buf_size-1); i >= 0; i--) {
        *buffer-- = (a & 1) + '0';
        a >>= 1;
    }
    for (size = buf_size-1;size >= 0;size--)
    {
	if (*(buffer+buf_size - size+1) == '1')
	    break;
    }
    return buffer;
}


void matrix_power(double* A,int N, int p, double* B)
{
    // p must >= 0
    // A B are complex matrices
    typedef struct{ double re; double im; } complex16;
    complex16 alpha;
    alpha.re = 1.0;
    alpha.im = 0.0;
    complex16 beta;
    beta.re = 0.0;
    beta.im = 0.0;
    if (p == 0)
    {
	# pragma omp parallel for
	for (int i=0; i < N; ++i)
	{
	    for (int j = 0; j< N; ++j)
	    {
		if (i==j)
		{
		    B[index(i,j,0,N)] = 1.0;
		    B[index(i,j,0,N)+1] = 0.0;
		}
		else
		{
		    B[index(i,j,0,N)] = 0.0;
		    B[index(i,j,0,N)+1] = 0.0;
		}
	    }
	}
    }
    else if (p <= 3)
    {
	double tmp[N*N*2];
	cblas_zcopy(N*N,A,1,B,1);
	for (int i =0; i< p-1;++i)
	{
	    cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,N,N,N,&alpha,B,N,A,N,&beta,tmp,N);
	    cblas_zcopy(N*N,tmp,1,B,1);
	}
    }
    else
    {
	char buffer[BUF_SIZE];
	buffer[BUF_SIZE - 1] = '\0';
	int size=0;
	int2bin(p, buffer, BUF_SIZE - 1,size);
	double* tmpptr,*tmp,*Z;
	Z = (double*) malloc (N*N*2*sizeof(double));
	tmp = (double*) malloc (N*N*2*sizeof(double));
	cblas_zcopy(N*N,A,1,Z,1);
	int q = 0;
	while (buffer[BUF_SIZE-2-q] == '0')
	{
	    cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,N,N,N,&alpha,Z,N,Z,N,&beta,tmp,N);
	    tmpptr = Z;
	    Z = tmp;
	    tmp = tmpptr;
	    q++;
	}
	cblas_zcopy(N*N,Z,1,B,1);	    
	for (int i=q+1;i<size;i++)
	{
	    cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,N,N,N,&alpha,Z,N,Z,N,&beta,tmp,N);
	    tmpptr = Z;
	    Z = tmp;
	    tmp = tmpptr;
	    if (buffer[BUF_SIZE-2-i]=='1')
	    {
		cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,N,N,N,&alpha,B,N,Z,N,&beta,tmp,N);
		cblas_zcopy(N*N,tmp,1,B,1);
	    }
	}
	free(Z);
	free(tmp);
    }
}

void matrix_vector_power(double* A,int N, int p, double* v_r, double* v)
{
    // power is 2^p
    // p must >= 2
    // A B are complex matrices
    typedef struct{ double re; double im; } complex16;
    complex16 alpha;
    alpha.re = 1.0;
    alpha.im = 0.0;
    complex16 beta;
    beta.re = 0.0;
    beta.im = 0.0;
    
    int k = (int)ceil(log2((double)N/log(2)));
    if (p<=k)
    {
	double* tmpptr,*tmp,*Z;
	Z = (double*) malloc (N*2*sizeof(double));
	tmp = (double*) malloc (N*2*sizeof(double));
	cblas_zcopy(N*2,v,1,Z,1);	    
	for (int i = 0; i < 2^p; i++)
	{
	    cblas_zgemv(CblasRowMajor,CblasNoTrans,N,N,&alpha,A,N,Z,1,&beta,tmp,1);
	    tmpptr = Z;
	    Z = tmp;
	    tmp = tmpptr;
	}
	cblas_zcopy(N*2,Z,1,v_r,1);	    
	free(Z);
	free(tmp);	
    }
    else
    {
	double* tmpptr,*tmp,*Z,*B;
	B = (double*) malloc (N*N*2*sizeof(double));
	Z = (double*) malloc (N*2*sizeof(double));
	tmp = (double*) malloc (N*2*sizeof(double));
	cblas_zcopy(N*2,v,1,Z,1);	    
	matrix_power(A,N,2^(p-k),B);
	for (int i = 0; i < 2^k; i++)
	{
	    cblas_zgemv(CblasRowMajor,CblasNoTrans,N,N,&alpha,A,N,Z,1,&beta,tmp,1);
	    tmpptr = Z;
	    Z = tmp;
	    tmp = tmpptr;
	}
	cblas_zcopy(N*2,Z,1,v_r,1);	    
	free(Z);
	free(tmp);	
	free(B);	
    }


}
void expm(double* A,int N,double t, double* B)
{
    
}

void pade3(double* A,double* U, double* V, int N)
{
    double b [] = {120., 60., 12., 1.};
}

void pade5(double* A,double* U, double* V, int N)
{
    double b [] = {30240., 15120., 3360., 420., 30., 1.};
}

void pade7(double* A,double* U, double* V, int N)
{
    double b [] = {17297280., 8648640., 1995840., 277200., 25200., 1512., 56., 1.};
}

void pade9(double* A,double* U, double* V, int N)
{
    double b [] = {17643225600., 8821612800., 2075673600., 302702400., 30270240.,2162160., 110880., 3960., 90., 1.};
}

void pade13(double* A,double* U, double* V, int N)
{
    double b [] = {64764752532480000., 32382376266240000., 7771770303897600., 1187353796428800., 129060195264000., 10559470521600., 670442572800.,33522128640., 1323241920., 40840800., 960960., 16380., 182., 1.};

}

*/
