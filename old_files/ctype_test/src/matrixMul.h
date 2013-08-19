#ifndef _MATRIXMUL_H_
#define _MATRIXMUL_H_
// Utilities and system includes
#include <stdio.h>
//define forward
/* extern "C" */
/* void deviceVerify(); */
extern "C"
void zdot(double* A,double* B,int N,double* result);
/* extern "C" */
/* void cdot(float* A,float* B,int N,float* result); */
#endif //#ifndef _MATRIXMUL_H_
