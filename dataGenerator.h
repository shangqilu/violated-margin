#ifndef DATAGENERATOR_H_INCLUDED
#define DATAGENERATOR_H_INCLUDED

#include "headers.h"

/*
*   generating a random number from rangeStart to rangeEnd
*/
double GenUniformRandom(double rangeStart, double rangeEnd);

/*
*   generating a random number from standard normal distribution
*/
double GenGaussianRandom();


/*
*   generating a Point with each dimensionality ranging from rangeStart to rangeEnd
*/
void GenRandomPoint(double rangeStart, double rangeEnd, Point &p);


/*
*   generating a Point in a d ball
*/
void GenRandomPointInCircle(Point &center, double radius, Point &p);

/*
*   generating a Hyperplane with each dimensionality ranging from rangeStart to rangeEnd
*/
void GenRandomHyperPlane(double rangeStart, double rangeEnd, HyperPlane &plane);


/*
*   generating a two dimension dataset
*	in the odd quadrant: label is 1 
*	in the even quadrant: label is -1
*/
void GenTwoDimensinoGridDataSet(char *filename, int totalNum);

/*
*   generating a Point in a ball and lies on the margin
*/
void GenMarginPoint(HyperPlane &plane, Point &center, double radius, double margin, Point &p);


/*
*   generating a dataset
*	all points are in a ball with given radius
*	and the optimal plane has the given margin
*/
void GenMarginDataSet(char *filename, double margin, double radius,
			int totalNum, int noiseNum);






#endif // DATAGENERATOR_H_INCLUDED
