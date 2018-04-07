#ifndef __PERCEPTRON_H__

#define __PERCEPTRON_H__
#include "headers.h"
#include <list>

/*
*   the max iterations in simple perceptron algorithm
*/
const int MaxIterations = 1000000;

/*
*   the simple perceptron algorithm for test
*/
void SimplePerceptron(PointSet &trainPoints, int dimension);

/*
*   compute a hyperplane with margin at least (1-epsilon)y_guess if y_guess < y_opt
*   R and y_min is the largest and smallest distance between a point 
*	from Points Set to origin in d+1 dimension
*/
bool BallMarginPerceptron(PointList &pointlist, HyperPlane &plane, int dimension, double y_guess, double R, double y_min);

/*
*   compute a separation hyperplane with margin at least (1-rho)*y_opt
*/
bool IncreMarginPerceptron(PointSet &trainPoints, HyperPlane &plane, int dimension, double rho);

#endif // PERCEPTRON_H
