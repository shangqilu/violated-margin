#ifndef __HEADERS_H__

#define __HEADERS_H__

#include <cstring>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
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

double Dot(double* w, double *x, int dimension);
void PrintHyperPlane(HyperPlane plane, int dimension);
double Distance(HyperPlane plane, Point pt, int dimension);
double Distance(Point pt1, Point pt2, int dimension);

#endif  //__HEADERS_H__
