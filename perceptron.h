#ifndef __PERCEPTRON_H__

#define __PERCEPTRON_H__
#include "headers.h"

/*
*   the max iterations in simple perceptron algorithm
*/
const int MaxIterations = 1000000;

/*
*   the simple perceptron algorithm for test
*/
void SimplePerceptron(PointSet &trainPoints, int dimension);

/*
*   compute a hyperplane with margin at least y_guess if y_guess < y_opt
*/
bool BallMarginPerceptron(PointList &pointlist, HyperPlane &plane, double y_guess);

/*
*   compute a separation hyperplane with margin at least (1-rho)*y_opt
*/
bool IncreMarginPerceptron(PointSet &trainPoints, HyperPlane &plane);

#endif // PERCEPTRON_H
