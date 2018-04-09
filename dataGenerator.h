#ifndef DATAGENERATOR_H_INCLUDED
#define DATAGENERATOR_H_INCLUDED

#include "headers.h"

/*
*   generating a random number from an interval
*/
double GenUniformRandom(double rangeStart, double rangeEnd);

/*
*   generating a random number from standard normal distribution
*/
double GenGaussianRandom();


void GenRandomPoint(double rangeStart, double rangeEnd, Point &p);

void GenRandomPointInCircle(Point &center, double radius, Point &p);

void GenRandomHyperPlane(double rangeStart, double rangeEnd, HyperPlane &plane);


void GenTwoDimensinoGridDataSet(char *filename, int totalNum);


void GenMarginPoint(HyperPlane &plane, Point &center, double radius, double margin, Point &p);

void GenMarginDataSet(char *filename, double margin, double radius,
			int totalNum, int noiseNum);






#endif // DATAGENERATOR_H_INCLUDED
