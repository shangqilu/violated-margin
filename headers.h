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
#include <list>

#define __DEBUG__
#ifdef __DEBUG__
#define DEBUG(format,...) printf("File: "__FILE__", Line: %05d: "format"\n", __LINE__, ##__VA_ARGS__)
#else
#define DEBUG(format,...)
#endif

#define PRINT_ERROR(format,...) printf("File: "__FILE__", Line: %05d: "format"\n", __LINE__, ##__VA_ARGS__)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

const double MAX_DOUBLE = 1e300;
const double ZERO = 0;
const double ZERO_ERROR = 1e-100;

/*
*   the structure of a hyperplane
*/
class HyperPlane
{
public:
	double *w;  //weights
	double b;   //bias
	int d;      //dimension
	HyperPlane(int d)
	{
		this->w = new double[d];
		for (int i = 0; i < d; i++) w[i] = 0;
		this->b = 0;
		this->d = d;
	}
	~HyperPlane()
	{
		delete []w;
	}
	
    void Clear()
    {
        for (int i = 0; i < this->d; i++)
            this->w[i] = 0;
        this->b = 0;
    }

private:
	HyperPlane(const HyperPlane &obj){}
	HyperPlane& operator = (const HyperPlane &obj){}

};

/*
*   the class of a point
*/
class Point{
public:
    double *x;
    int y;
    int dim;
    Point(int d)
    {
        x = new double[d];
        for (int i = 0; i < d; i++)
            x[i] = 0;
        y = 0;
        this->dim = d;
    }
    ~Point()
    {
        if (x != NULL)
        {
            delete []x;
            x = NULL;
        }
    }

    Point(const Point& obj) 
    {
        x = new double[obj.dim];
        for (int i = 0; i < obj.dim; i++)
        {
            x[i] = obj.x[i];
        }
        y = obj.y;
        dim = obj.dim;
    }
    Point& operator = (const Point & obj) 
	{
        if (this == &obj) {
            return *this;
        }
        else {
            delete[]x;
            x = new double[obj.dim];
            for (int i = 0; i < obj.dim; i++)
            {
                x[i] = obj.x[i];
            }
            y = obj.y;
            dim = obj.dim;
            return *this;
        }
    }
private:
};

/*
*   Points Set
*/
typedef vector<Point> PointSet;
/*
*   Points List
*/
typedef list<Point> PointList;
/*
*	Point Index
*/
typedef vector<int> PointIndex;

/*
*   Load data from two files, a data file and a label file
*/
PointSet LoadData(char* filename, char* label_filename, int dimension);

/*
*   load data from single file in the format of libSVM
*/
PointSet LoadDataLibSVMFormat(char* filename, int dimension);

/*
*   return a Point Set by copying the data set
*/
PointSet CopyPoints(PointSet &points, int dimension);

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
/*
*   compute the vector of p1-p2
*/
Point PointMinus(Point &pt1, Point &pt2, int dimension);

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
bool MinimumViolatedDistance(PointSet &points, HyperPlane &plane, double &min_dis, int k, int &reak_k);

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


/*
*   Print Error
*/
void PrintError(string msg);

#endif  //__HEADERS_H__
