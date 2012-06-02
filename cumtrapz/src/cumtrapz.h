#ifndef _CUMTRAPZ_H_
#define _CUMTRAPZ_H_
// Utilities and system includes
#include <stdio.h>
#include "math.h"
#include "mkl.h"
#include "mkl_cblas.h"
#define BUF_SIZE 50 		/* used in int2bin */
//define forward
/* extern "C" */
/* void deviceVerify(); */
int inline index(int i,int j,int k);
extern "C"
void addorder(double* T,double* D,double* env, double* A,int N, int num,double dt);
extern "C"
void cumtrapz(double* A,double dt,int N,int num);
extern "C"
void matrix_power(double* A,int N, int p, double* B);
char *int2bin(int a, char *buffer, int buf_size, int& size);
extern "C"
void expm(double* A,int N,double t,double* B);
void pade3(double* A,double* U, double* V, int N);
void pade5(double* A,double* U, double* V, int N);
void pade7(double* A,double* U, double* V, int N);
void pade9(double* A,double* U, double* V, int N);
void pade13(double* A,double* U, double* V, int N);
/* extern "C" */
/* void cdot(float* A,float* B,int N,float* result); */
#endif //#ifndef _MATRIXMUL_H_
