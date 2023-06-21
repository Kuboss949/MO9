//
// Created by creaz on 09.05.2023.
//

#include "thomas.h"

void matrixProcedure(const double * lowDiagonal, double * diagonal, const double * upDiagonal, int n){
    for(int i = 1; i<n; i++){
        diagonal[i] = diagonal[i] - lowDiagonal[i-1]*(1.0/diagonal[i-1])*upDiagonal[i-1];
    }
}

void vectorProcedure(double *b, const double * lowDiagonal, const double * diagonal, const double * upDiagonal, int n){
    for(int i = 1; i<n; i++){
        b[i] = b[i] - lowDiagonal[i-1]*(1.0/diagonal[i-1])*b[i-1];
    }
    b[n-1] = b[n-1]*1.0/diagonal[n-1];
    for(int i = n-2; i>=0; i--){
        b[i] = (1.0/diagonal[i])*(b[i]-upDiagonal[i]*b[i+1]);
    }
}
