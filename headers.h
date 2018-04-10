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

/*
*   sometimes there can a trouble in simplex algorithm
*	if ZERO_ERROR is very small, may be there would be no solution
*	if ZERO_ERROR is very large, may be there would be a cycling
*	
*
*/
const double ZERO_ERROR = 1e-15;

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
    Point(int d, double *x)
    {
		this->x = x;
        for (int i = 0; i < d; i++)
            x[i] = 0;
        y = 0;
    }
private:
	//Point(const Point& obj){}
	//Point& operator = (const Point & obj) {}
};

/*
*   Points Set
*/
typedef vector<Point*> PointSet;
/*
*   Points List
*/
typedef list<Point*> PointList;
/*
*	Point Index
*/
typedef vector<int> PointIndex;





/*
*   Print Error
*/
void PrintError(string msg);

#endif  //__HEADERS_H__
