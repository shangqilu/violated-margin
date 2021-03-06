#ifndef TOOLS_H_INCLUDED
#define TOOLS_H_INCLUDED


#include "headers.h"



/*
*   copy the content of plane2 to plane1
*/
void CopyHyperPlane(HyperPlane &plane1, HyperPlane &plane2);

/*
*   compute the dot product of two array
*/
double Dot(double* w, double *x, int dimension);

/*
*   compute the dot product of two points
*/
double Dot(Point &pt1, Point &pt2, int dimension);

/*
*   compute the dot product of two vectors
*/
double Dot(vector<double> x, vector<double> y, int dimension);

void PrintHyperPlane(HyperPlane &plane);
/*
*   compute the distance from a point to a plane
*/
double Distance(HyperPlane &plane, Point &pt, int dimension);
/*
*   compute the distance between two points
*/
double Distance(Point &pt1, Point &pt2, int dimension);


/*
*   return whether points can be separated correctly by a hyperplane
*       if return true min_dis record the minimal distance from points to the plane
*       else return false
*/
bool MinimumSeparableDistance(PointSet &points, HyperPlane &plane, double &min_dis);

/*
*   return whether sub set points can be separated by a hyperplane without violating k points
*       if return true min_dis records the minimal distance from points to the plane
*       else return false real_k records the number of points the plane violates
*/
bool MinimumViolatedDistance(PointSet &points, HyperPlane &plane, double &min_dis, int k, int &real_k);

/*
*   return whether points can be separated by a hyperplane without violating k points
*       if return true min_dis records the minimal distance from points to the plane
*       else return false real_k records the number of points the plane violates
*/
bool MinimumSubsetViolatedDistance(PointSet &points, PointIndex &index,
	HyperPlane &plane, double &min_dis, int k, int &real_k);

/*
*   compute the solution of a equation set
*   if return true array x contains the solution
*/
bool GaussianEquation(double** A, double* b, double* x, int n);
/*
*   compute the inverse of Matrix A
*   if return true the inverse Matrix is stored in B
*/
bool GaussianInverseMatrix(double** A, double**B, int n);


void PrintMatrix(double **T, int dimension);

/*
*   compute Matrix A * Matrix B
*   the result is stored in C
*/
void MatrixMultiply(double **A, double **B, double **C, int dimension);

/*
*   compute the transformation of a point set
*   T is the transformation Matrix which has size (d+1)*(d+1)
*   d is the dimension of points
*/
void TransformingPoints(PointSet &points, double **T, int point_dimension);

void PrintPoints(PointSet &points, int dimension);

/*
*   compute the largest prime which is no more greater than m
*/
int LargestPrimeBelow(int m);
/*
*   determine whether n is prime
*/
bool is_prime(int n);
/*
*   compute the factorial of d
*/
long long Factorial(int d);

/*
*   compute the Log Factorial of d
*/
double LogFactorial(int d);

#endif // TOOLS_H_INCLUDED
