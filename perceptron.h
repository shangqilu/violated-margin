#ifndef __PERCEPTRON_H__

#define __PERCEPTRON_H__
#include "headers.h"


const int MaxIterations = 1000000;
void SimplePerceptron(PointSet trainPoints, int dimension);
bool MarginPerceptrion(PointSet trainPoints, HyperPlane &plane, int dimension, double y_guess, double R, double epsilon);
HyperPlane IncreMarginPerceptrion(PointSet trainPoints, int dimension, double rho);

#endif // PERCEPTRON_H
