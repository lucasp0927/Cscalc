#ifndef _CUMTRAPZ_H_
#define _CUMTRAPZ_H_
// Utilities and system includes
#include <stdio.h>
#include "mkl.h"
#include "mkl_cblas.h"
//define forward
/* extern "C" */
/* void deviceVerify(); */
extern "C"
void addorder(double* T,double* D,double* env, double* A,int N, int num,double dt);
extern "C"
void cumtrapz(double* A,double dt,int N,int num);
int inline index(int i,int j,int k);
/* extern "C" */
/* void cdot(float* A,float* B,int N,float* result); */
#endif //#ifndef _MATRIXMUL_H_
