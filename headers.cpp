#include "headers.h"

using namespace std;


double Dot(double* w, double *x, int dimension)
{
    double sum = 0;
    for (int i = 0; i < dimension; i++)
    {
        sum += (*w) * (*x);
        w++;
        x++;
    }
    return sum;
}

double Dot(Point pt1, Point pt2, int dimension)
{
    double sum = 0;
    for (int i = 0; i < dimension; i++)
    {
        sum += pt1.x[i]*pt2.x[i];
    }
    return sum;
}

Point PointMinus(Point pt1, Point pt2, int dimension)
{
    Point ans = Point(dimension);
    for (in i = 0; i < dimension; i++)
    {
        ans.x[i] = pt1.x[i] - pt2.x[i];
    }
    return ans;
}

void PrintHyperPlane(HyperPlane plane, int dimension)
{
    cout << "w: (";
    for (int i = 0; i < dimension; i++)
    {
        cout << plane.w[i];
        if (i != dimension - 1)
        {
            cout << ", ";
        }
    }
    cout << ") ";
    cout << "b: " << plane.b << endl;
}

double Distance(HyperPlane plane, Point pt, int dimension)
{
    double dis = 0, denominator = 0;
    for (int i = 0; i < dimension; i++)
    {
        dis += plane.w[i]*pt.x[i];
        denominator += plane.w[i]*plane.w[i];
    }
    dis += plane.b;
    return fabs(dis)/sqrt(denominator);
}

double Distance(Point pt1, Point pt2, int dimension)
{
    double dis = 0;
    for(int i = 0; i < dimension; i++)
    {
        dis += (pt1.x[i] - pt2.x[i]) * (pt1.x[i] - pt2.x[i]);
    }
    return sqrt(dis);
}



bool Gaussian(double** A, double* b, double* x, int n)
{


    double** M = new double*[n];
    for (int i=0; i<n; i++)
    {
        M[i] = new double[n+1];
        for (int j=0; j<n; j++)
        {
            M[i][j] = A[i][j];
        }
        M[i][n] = b[i];
    }
    for (int k=0; k<n; k++)
    {
        //
        double colMax = fabs(M[k][k]);
        int maxLineIndex = k;
        for(int i=k+1; i<n; i++)
        {
            //find the largest element in current column
            if(fabs(M[i][k]) > colMax)
            {
                colMax = fabs(M[i][k]);
                maxLineIndex = i;
            }
        }
        if(colMax < EPSILONG)
        {
            //singular matrix
            for (int i=0; i<n; i++)
            {
                delete M[i];
            }
            delete M;
            return false;
        }
        double temp;
        //swap line k and maxLineIndex
        for (int m=0; m<n+1; m++)
        {
            temp = M[k][m];
            M[k][m] = M[maxLineIndex][m];
            M[maxLineIndex][m] = temp;
        }

        //
        for(int i=k+1; i<n; i++)
        {
            for (int j=k+1; j<n+1; j++)
            {
                M[i][j] = M[k][k]*M[i][j]/M[i][k] - M[k][j];
            }
        }
    }
    //calculate the solution
    for (int i=n-1; i>=0; i--)
    {
        x[i] = M[i][n];
        for (int j=i+1; j<n; j++)
        {
            x[i] -= M[i][j]*x[j];
        }
        x[i] /= M[i][i];
    }
    for (int i=0; i<n; i++)
    {
        delete M[i];
    }
    delete M;
    return true;
}
