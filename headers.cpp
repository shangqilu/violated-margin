#include "headers.h"

using namespace std;


double Dot(double* w, double *x, int dimension)
{
    double sum = 0;
    for (int i = 0; i < dimension; i++) {
        sum += (*w) * (*x);
        w++;
        x++;
    }
    return sum;
}

void PrintHyperPlane(HyperPlane plane, int dimension)
{
    cout << "w: (";
    for (int i = 0; i < dimension; i++) {
        cout << plane.w[i];
        if (i != dimension - 1){
            cout << ", ";
        }
    }
    cout << ")";
    cout << "b: " << plane.b << endl;
}

double Distance(HyperPlane plane, Point pt, int dimension)
{
    double dis = 0, denominator = 0;
    for (int i = 0; i < dimension; i++) {
        dis += plane.w[i]*pt.x[i];
        denominator += plane.w[i]*plane.w[i];
    }
    dis += plane.b;
    return fabs(dis)/sqrt(denominator);
}

double Distance(Point pt1, Point pt2, int dimension)
{
    double dis = 0;
    for(int i = 0; i < dimension; i++) {
        dis += (pt1.x[i] - pt2.x[i]) * (pt1.x[i] - pt2.x[i]);
    }
    return sqrt(dis);
}
