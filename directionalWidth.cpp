#include "directionalWidth.h"
#include <ctime>


using namespace std;

void TwoApproxiDiameter(PointSet points, int dimension, int &s, int &t)
{
    //compute the 2-approximation diameter
    srand(time(NULL));
    while(1)
    {
        int s = rand()% n;
        double max_dis = MAX_DOUBLE;
        double t = -1;
        for (int i = 0; i < n; i++)
        {
            if (i != s)
            {
                double cur_dis = Distance(points[s], points[i], dimension);
                if (cur_dis > max_dis)
                {
                    t = i;
                    max_dis = cur_dis;
                }
            }
        }
        int flag = 0;
        Point direction = PointMinus(points[s], points[t], dimension);
        Point axis = Point(dimension);
        axis.x[dimension-1] = 1;
        //the hyper plane s and t defined is not parallel to last dimension coordinate axis
        if (fabs(Dot(direction, axis)) > ZERO)
        {
            break;
        }
    }
}


void MinimumBoudingBox(PointSet points, int dimension)
{
    int n = points.size();
    if (n < 2)
    {
        puts("there are no more than 2 points!");
        return;
    }
    int s, t;
    TwoApproxiDiameter(points, dimension, s, t);
    //all points is above point t
    Point d_axis = PointMinus(points[s], points[t], dimension);
    //find those d-1 basis


    double **CurT = new double*[dimension+1];
    //store those d basis vectors
    for (int i = 0; i < dimension+1; i++)
    {
        CurT = new double[dimension+1];
        for (int j = 0; j < dimension+1; j++)
        {
            CurT[i][j] = 0;
        }
    }
    for (int i = 0; i < dimension; i++)
    {
        CurT[i][0] = d_axis.x[i];
    }


    double **A = new double*[dimension-1];
    for (int i = 0; i < dimension-1; i++)
    {
        A[i] = new double[dimension-1];
    }
    double *b = new double[dimension-1];
    double *x = new double[dimension-1];
    for (int number_variables = 1; number_variables <= dimension-1; number_variables++)
    {
        // there are number_variables variables
        //set the previous dimension-i variables to 1
        //form the equation set
        for (int i = 0; i < number_variables; i++)
        {
            double sum = 0;
            for (int j = 0; j < dimension - number_variables; j++)
            {
                sum += CurT[j][i]*(1-points[t].x[j]);
            }
            for (int j = dimension - number_variables; j < dimension; j++)
            {
                A[i][j- (dimension-number_variables)] = CurT[j][i];
                sum -= CurT[j][i] * points[t].x[j];
            }
            b[i] = -sum;
        }
        bool solve = Gaussian(A, b, x, number_variables);
        if (!solve)
        {
            puts("the equation set is insouble!");
            break;
        }

        for (int j = 0; j < dimension - number_variables; j++)
        {
            CurT[j][number_variables] = 1;
        }
        for (int j = dimension - number_variables; j < dimension; j++)
        {
            cur[j][number_variables] = x[j-(dimension-number_variables)]
        }
    }




}
