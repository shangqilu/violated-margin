#include "perceptron.h"
#include <iostream>

using namespace std;



void SimplePerceptron(PointSet trainPoints, int dimension)
{
    int n = trainPoints.size();
    HyperPlane plane = HyperPlane(dimension);
    int iter_cnt = 0;
    for (iter_cnt = 0; iter_cnt < MaxIterations; iter_cnt++)
    {
        printf("Iteration %d\n", iter_cnt);
        int i;
        for(i = 0; i < n; i++)
        {
            Point cur_pt = trainPoints[i];
            if(cur_pt.y * (Dot(plane.w, cur_pt.x, dimension) + plane.b)   <= 0)
            {
                for (int j = 0; j < dimension; j++)
                {
                    plane.w[j] += cur_pt.y*cur_pt.x[j];
                }
                plane.b += cur_pt.y;
                break;
            }
        }
        if(i == n)
            break;
    }
    if (iter_cnt == MaxIterations) {
        puts("It is non-linearly separable!");
    }
    else PrintHyperPlane(plane, dimension);
}

bool MarginPerceptron(PointSet trainPoints, HyperPlane &plane, int dimension, double y_guess, double R, double epsilon)
{
    int n = trainPoints.size();
    int iter_cnt = 0;
    double maxIters = (2+2*epsilon)*R*R/(epsilon*epsilon*y_guess*y_guess);
    printf("maxIters: %lf\n", maxIters);
    cout << "R: " << R << " y_guess: " << y_guess << endl;
    for (iter_cnt = 0; iter_cnt < maxIters; iter_cnt++)
    {
        /*
        if (iter_cnt % 500 == 0)
        {
            printf("Iteration %d\n", iter_cnt);
        }
        */
        int i;
        for(i = 0; i < n; i++)
        {
            Point cur_pt = trainPoints[i];
            if(Distance(plane, cur_pt, dimension) <= (1-epsilon)*y_guess
               || cur_pt.y * (Dot(plane.w, cur_pt.x, dimension) + plane.b)   <= 0)
            {
                for (int j = 0; j < dimension; j++)
                {
                    plane.w[j] += cur_pt.y*cur_pt.x[j];
                }
                plane.b += cur_pt.y;
                break;
            }
        }
        if(i == n)
            break;
    }
    if (iter_cnt >= maxIters) {
        puts("It is non-linearly separable based on y_guess");
        return false;
    }
    else {
        return true;
    }
}

bool IncreMarginPerceptron(PointSet trainPoints, HyperPlane &plane, int dimension, double rho)
{

    //find the max distance between the points and the origin
    double maxDistantce = 0, minDistance = MAX_DOUBLE;
    int n = trainPoints.size();

    double epsilon = 1 - sqrt(1 - rho);   //parameter in margin perceptron
    Point origin = Point(dimension);
    for (int i = 0; i < n; i++) {
        double curDis = Distance(trainPoints[i], origin, dimension);
        if (maxDistantce < curDis) maxDistantce = curDis;
        if (minDistance > curDis) minDistance = curDis;
    }
    cout << "Min: " << minDistance << "  Max: " << maxDistantce << endl;
    double y_guess = maxDistantce;
    while(true) {
        plane.Clear();
        bool re = MarginPerceptron(trainPoints, plane, dimension, y_guess, maxDistantce, epsilon);
        if (re) {
            break;
        } else{
            y_guess = (1-epsilon)*y_guess;
        }
        if(y_guess < 0.1 ) {////need to consider carefully
            puts("It is non-linear separable!");
            return false;
        }
    }
    return true;
}
