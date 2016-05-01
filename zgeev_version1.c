/*******************************************************************************
 *  Copyright (C) 2009-2015 Intel Corporation. All Rights Reserved.
 *  The information and material ("Material") provided below is owned by Intel
 *  Corporation or its suppliers or licensors, and title to such Material remains
 *  with Intel Corporation or its suppliers or licensors. The Material contains
 *  proprietary information of Intel or its suppliers and licensors. The Material
 *  is protected by worldwide copyright laws and treaty provisions. No part of
 *  the Material may be copied, reproduced, published, uploaded, posted,
 *  transmitted, or distributed in any way without Intel's prior express written
 *  permission. No license under any patent, copyright or other intellectual
 *  property rights in the Material is granted to or conferred upon you, either
 *  expressly, by implication, inducement, estoppel or otherwise. Any license
 *  under such intellectual property rights must be express and approved by Intel
 *  in writing.
 *
 ********************************************************************************
 */
/*
 ZGEEV Example.
 ==============
 
 Program computes the eigenvalues and left and right eigenvectors of a general
 rectangular matrix A:
 
 ( -3.84,  2.25) ( -8.94, -4.75) (  8.95, -6.53) ( -9.87,  4.82)
 ( -0.66,  0.83) ( -4.40, -3.82) ( -3.50, -4.26) ( -3.15,  7.36)
 ( -3.99, -4.73) ( -5.88, -6.60) ( -3.36, -0.40) ( -0.75,  5.23)
 (  7.74,  4.18) (  3.66, -7.53) (  2.58,  3.60) (  4.59,  5.41)
 
 Description.
 ============
 
 The routine computes for an n-by-n complex nonsymmetric matrix A, the
 eigenvalues and, optionally, the left and/or right eigenvectors. The right
 eigenvector v(j) of A satisfies
 
 A*v(j)= lambda(j)*v(j)
 
 where lambda(j) is its eigenvalue. The left eigenvector u(j) of A satisfies
 
 u(j)H*A = lambda(j)*u(j)H
 
 where u(j)H denotes the conjugate transpose of u(j). The computed
 eigenvectors are normalized to have Euclidean norm equal to 1 and
 largest component real.
 
*/


/* Comple this code as gcc -O zgeev_version1.c -o ssy -llapack -lblas */



#include <stdlib.h>
#include <stdio.h>

/* Complex datatype */
struct _dcomplex { double re, im; };
typedef struct _dcomplex dcomplex;

/* ZGEEV prototype */
extern void zgeev( char* jobvl, char* jobvr, int* n, dcomplex* a,
                  int* lda, dcomplex* w, dcomplex* vl, int* ldvl, dcomplex* vr, int* ldvr,
                  dcomplex* work, int* lwork, double* rwork, int* info );
/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, dcomplex* a, int lda );

/* Parameters */
#define N 4
#define LDA N
#define LDVL N
#define LDVR N

/* Main program */
int main(void) {
    /* Locals */
    int n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR, info, lwork;
    dcomplex wkopt;
    dcomplex* work;
    /* Local arrays */
    /* rwork dimension should be at least 2*n */
    double rwork[2*N];
    dcomplex w[N], vl[LDVL*N], vr[LDVR*N];
    dcomplex a[LDA*N] = {
        {-3.84,  2.25}, {-0.66,  0.83}, {-3.99, -4.73}, { 7.74,  4.18},
        {-8.94, -4.75}, {-4.40, -3.82}, {-5.88, -6.60}, { 3.66, -7.53},
        { 8.95, -6.53}, {-3.50, -4.26}, {-3.36, -0.40}, { 2.58,  3.60},
        {-9.87,  4.82}, {-3.15,  7.36}, {-0.75,  5.23}, { 4.59,  5.41}
    };
    /* Executable statements */
    printf( " ZGEEV Example Program Results\n" );
    /* Query and allocate the optimal workspace */
    lwork = -1;
    zgeev( "Vectors", "Vectors", &n, a, &lda, w, vl, &ldvl, vr, &ldvr,
          &wkopt, &lwork, rwork, &info );
    lwork = (int)wkopt.re;
    work = (dcomplex*)malloc( lwork*sizeof(dcomplex) );
    /* Solve eigenproblem */
    zgeev( "Vectors", "Vectors", &n, a, &lda, w, vl, &ldvl, vr, &ldvr,
          work, &lwork, rwork, &info );
    /* Check for convergence */
    if( info > 0 ) {
        printf( "The algorithm failed to compute eigenvalues.\n" );
        
    }
    /* Print eigenvalues */
    print_matrix( "Eigenvalues", 1, n, w, 1 );
    /* Print left eigenvectors */
    print_matrix( "Left eigenvectors", n, n, vl, ldvl );
    /* Print right eigenvectors */
    print_matrix( "Right eigenvectors", n, n, vr, ldvr );
    /* Free workspace */
    free( (void*)work );
    return 0;
} /* End of ZGEEV Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, dcomplex* a, int lda ) {
    int i, j;
    printf( "\n %s\n", desc );
    for( i = 0; i < m; i++ ) {
        for( j = 0; j < n; j++ )
        printf( " (%6.2f,%6.2f)", a[i+j*lda].re, a[i+j*lda].im );
        printf( "\n" );
    }
}