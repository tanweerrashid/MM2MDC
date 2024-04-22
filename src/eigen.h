#ifndef EIGEN_H
#define EIGEN_H

#define TOLERANCE 0.001f

//#include "GeoCommon.h"

/**
 * Numerical functions for computing minimizers of a least-squares system
 * of equations.
 *
 * @author Scott Schaefer
 */


/**
 * Uses a jacobi method to return the eigenvectors and eigenvalues 
 * of a 3x3 symmetric matrix.  Note: "a" will be destroyed in this
 * process.  "d" will contain the eigenvalues sorted in order of 
 * decreasing modulus and v will contain the corresponding eigenvectors.
 *
 * @param a the 3x3 symmetric matrix to calculate the eigensystem for
 * @param d the variable to hold the eigenvalues
 * @param v the variables to hold the eigenvectors
 */
void jacobi ( float a[][3], float d[], float v[][3] );

/**
 * Inverts a 3x3 symmetric matrix by computing the pseudo-inverse.
 *
 * @param mat the matrix to invert
 * @param midpoint the point to minimize towards
 * @param rvalue the variable to store the pseudo-inverse in
 * @param w the place to store the inverse of the eigenvalues
 * @param u the place to store the eigenvectors
 */
void matInverse ( float mat[][3], float midpoint[], float rvalue[][3], float w[], float u[][3] );

/**
 * Calculates the L2 norm of the residual (the error)
 * (Transpose[A].A).x = Transpose[A].B
 * 
 * @param a the matrix Transpose[A].A
 * @param b the matrix Transpose[A].B
 * @param btb the value Transpose[B].B
 * @param point the minimizer found
 *
 * @return the error of the minimizer
 */
float calcError ( float mat[], float point[] );

/**
 * Calculates the minimizer of the given system and returns its error.
 *
 * @param midpoint the point to minimize towards
 * @param rvalue the place to store the minimizer
 * @param box the volume bounding the voxel this QEF is for
 *
 * @return the error of the minimizer
 */
float calcPoint ( float midpoint[], float rvalue[], float *mat );

void qr ( float eqs[][4], int num, float *rvalue );
void qr ( float *mat, float *eqn );
void qr ( float *mat1, float *mat2, float *rvalue );
void qr ( float eqs[][4], int num = 4 );

int estimateRank ( float *a );

#endif