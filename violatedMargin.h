#ifndef VIOLATEDMARGIN_H_INCLUDED
#define VIOLATEDMARGIN_H_INCLUDED

#include "headers.h"
#include "perceptron.h"
#include "simplex.h"
#include "directionalWidth.h"

using namespace std;


/*
*   Sampling the point set with probability p
*   return the sample set
*/
PointSet Sampling(PointSet &points, int dimension, double p);

/*
*   choose a margin classification algorithm
*		method = 0: perceptron algorithm
*		method = 1: simplex algorithm
*		method = 2: directional width algorithm
*/
bool MarginClasification(PointSet &points, HyperPlane &plane, int dimension, double rho,
	int method, PointSet &coresetDirections, PointSet &classifyDirections);

/*
*	compute a hyperplane violated no more than (1+epsilon)*k points
*	and the margin is no less than (1- rho)optimal_margin
*	with high successful probability no less than 1 - delta
*/
bool ApproximateViolatedMargin(PointSet &points, HyperPlane &optimal_plane, int dimension,
                                int k, double epsilon, double rho, double delta, int method);

/*
*	compute a hyperplane violated no more than k points
*	and the margin is no less than (1- rho)optimal_margin
*	with high successful probability no less than 1 - delta
*/
bool ViolatedMargin(PointSet &points, HyperPlane &optimal_plane, int dimension,
                     int k, double rho, double delta, int method);


#endif // VIOLATEDMARGIN_H_INCLUDED
