#ifndef MATRIX_OPERATIONS_H 
#define MATRIX_OPERATIONS_H

#include <stdio.h>
#include "arrayFunctions.h"
#include <gmp.h>

typedef struct sMtx {
    int     dim_x;
    int     dim_y;
    mpz_t     *mtx;
} *Matrix, sMatrix;


Matrix calculateGaussianMatrix(VectorSpace * v, int rowSize, mpz_t pmenusone);
void MtxNormalizeRow( Matrix m, int rix, int lead, mpz_t pmenusone);
void MtxToReducedREForm(Matrix m, mpz_t pmenusone);
void MtxSwapRows( Matrix m, int rix1, int rix2);
void MtxMulAndAddRows(Matrix m, int ixrdest, int ixrsrc, mpz_t mplr, mpz_t pmenusone);
void MtxDisplay( Matrix m );
Matrix InitMatrix( int x_dim, int y_dim, VectorSpace * v);
void MtxSetRow(Matrix m, int irow, int *v);
Matrix NewMatrix( int x_dim, int y_dim );
int ownMatrixAlgorithm(Matrix m, mpz_t pmenusone);

#endif


