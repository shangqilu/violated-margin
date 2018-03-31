#ifndef VIOLATEDMARGIN_H_INCLUDED
#define VIOLATEDMARGIN_H_INCLUDED

#include "headers.h"
#include "perceptron.h"
#include "simplex.h"
#include "directionalWidth.h"

using namespace std;


//sampling with probability p
//copy those points into a new Point Set
PointSet Sampling(PointSet points, int dimension, double p);
bool MarginClasification(PointSet points, HyperPlane &plane, int dimension, double rho, int method);

bool ApproximateViolatedMargin(PointSet points, HyperPlane &optimal_plane, int dimension,
                                int k, double epsilon, double rho, double delta, int method);

bool ViolatedMargin(PointSet points, HyperPlane &optimal_plane, int dimension,
                     int k, double rho, double delta, int method);


#endif // VIOLATEDMARGIN_H_INCLUDED
