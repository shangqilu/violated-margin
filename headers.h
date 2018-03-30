#ifndef __HEADERS_H__

#define __HEADERS_H__

#include <cstring>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>


#define __DEBUG__
#ifdef __DEBUG__
#define DEBUG(format,...) printf("File: "__FILE__", Line: %05d: "format"\n", __LINE__, ##__VA_ARGS__)
#else
#define DEBUG(format,...)
#endif


using namespace std;

#define MAX_DOUBLE 1e9
#define ZERO  1e-6
struct HyperPlane
{
    double *w;
    double b;
    int d;
    HyperPlane(int d)
    {
        this->w = new double[d];
        for (int i = 0; i < d; i++)
            this->w[i] = 0;
        this->b = 0;
        this->d = d;
    }
    void Clear()
    {
        for (int i = 0; i < this->d; i++)
            this->w[i] = 0;
        this->b = 0;
    }
};

struct Point
{
    double *x;
    int y;
    Point(int d)
    {
        this->x = new double[d];
        for (int i = 0; i < d; i++)
            this->x[i] = 0;
        this->y = 0;
    }
    Point(int d, double *x, int y)
    {
        this->x = new double[d];
        for (int i = 0; i < d; i++)
            this->x[i] = x[i];
        this->y = y;
    }

};

typedef vector<Point> PointSet;

PointSet LoadData(char* filename, char* label_filename, int dimension);


//load data from single file in the format of libSVM
PointSet LoadDataLibSVMFormat(char* filename, int dimension);

PointSet CopyPoints(PointSet points, int dimension);
void CopyHyperPlane(HyperPlane &plane1, HyperPlane &plane2);

double Dot(double* w, double *x, int dimension);

double Dot(Point pt1, Point pt2, int dimension);

Point PointMinus(Point pt1, Point pt2, int dimension);

void PrintHyperPlane(HyperPlane plane, int dimension);

double Distance(HyperPlane plane, Point pt, int dimension);

double Distance(Point pt1, Point pt2, int dimension);
double MinimumSeparableDistance(PointSet points, HyperPlane plane);

bool GaussianEquation(double** A, double* b, double* x, int n);
bool GaussianInverseMatrix(double** A, double**B, int n);


void PrintMatrix(double **T, int dimension);
void MatrixMultiply(double **A, double **B, double **C, int dimension);
void TransformingPoints(PointSet &points, double **T, int point_dimension);

void PrintPoints(PointSet points, int dimension);


int LargestPrimeBelow(int m);
bool is_prime(int n);
long long Factorial(int d);
double LogFactorial(int d);

#endif  //__HEADERS_H__
