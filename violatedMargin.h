#ifndef VIOLATEDMARGIN_H_INCLUDED
#define VIOLATEDMARGIN_H_INCLUDED

#include "headers.h"


/*
*	compute a separation hyperplane violating no more than (1+epsilon)*k points
*	and the margin is no less than (1- rho)*optimal_margin
*/
bool ApproximationViolatedMargin(PointSet &points, HyperPlane &optimal_plane);
/*
*	compute a separation hyperplane violating no more than k points
*	and the margin is no less than (1- rho)optimal_margin
*/
bool ViolatedMargin(PointSet &points, PointIndex &subSetIndex, HyperPlane &optimalPlane, int k_prime);



/*
*   Sampling the original point set with index starting from 1 to n
*	with probability p
*   return the sampled set storing indexes
*/
PointIndex Sampling(int n, double p);
/*
*   Sampling a sub point set storing indexes of original set with probability p
*   return the sample set storing indexes
*/
PointIndex Sampling(PointIndex &points, double p);


/*
*   choose a margin classification algorithm
*		method = 0: perceptron algorithm
*		method = 1: simplex algorithm
*		methdo = 2: directional width algorithm
*/
bool MarginClasification(PointSet &points, HyperPlane &plane, PointSet &coresetDirections, PointSet &classifyDirections);
#endif // VIOLATEDMARGIN_H_INCLUDED
