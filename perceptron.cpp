#include "perceptron.h"
#include <iostream>

using namespace std;



void SimplePerceptron(PointSet &trainPoints, int dimension)
{
    int n = trainPoints.size();
    HyperPlane plane(dimension);
    int iter_cnt = 0;
    for (iter_cnt = 0; iter_cnt < MaxIterations; iter_cnt++)
    {
        printf("Iteration %d\n", iter_cnt);
        int i;
        for(i = 0; i < n; i++)
        {
            Point cur_pt = trainPoints[i];
            if(cur_pt.y * (Dot(plane.w, cur_pt.x, dimension) + plane.b) <= 0)
            {
                for (int j = 0; j < dimension; j++)
                {
                    plane.w[j] += cur_pt.y*cur_pt.x[j];
                }
                plane.b += cur_pt.y;
				printf("break %d %d\n", i, cur_pt.y);
                break;
            }
        }
        if(i == n)
            break;
    }
    if (iter_cnt == MaxIterations) {
        puts("It is non-linearly separable!");
    }
    else PrintHyperPlane(plane);
}


bool BallMarginPerceptron(PointSet &trainPoints, HyperPlane &plane, int dimension, double y_guess, double R, double y_min)
{
	int n = trainPoints.size();
	int iter_cnt = 0;
	int maxIters = (int)((R*R)/(y_min*y_min) + 0.5);
	maxIters = 10000;
	printf("maxIters: %d\n", maxIters);
	cout << "R: " << R << " y_guess: " << y_guess << " y_min" << y_min << endl;
	for (iter_cnt = 0; iter_cnt < maxIters; iter_cnt++)
	{
		int i = 0;
		for (i = 0; i < n; i++)
		{
			Point cur_pt = trainPoints[i];
			//check whether (d+1) dimension cur plane passes thorough 
			//d dimension Ball(points[i], y_guess);
			if (cur_pt.y * (Dot(plane.w, cur_pt.x, dimension) + plane.b) < 0)
			{
				double denominator = 0;
				for (int j = 0; j < dimension; j++) {
					denominator += plane.w[j];
				}
				//W_t+1 = W_t+ y*(x - y*y_guess*W_t/||W_t||)
				for (int j = 0; j < dimension; j++)
				{
					plane.w[j] += cur_pt.y*(cur_pt.x[j] - 
						cur_pt.y * y_guess * plane.w[j]/sqrt(denominator));
				}
				plane.b += cur_pt.y;
				printf("break %d\n", i);
				break;
			}
		}
		if (i == n)
		{
			break;
		}
	}
	if (iter_cnt >= maxIters) {
		puts("It is non-linearly separable based on y_guess");
		return false;
	}
	else {
		PrintHyperPlane(plane);
		return true;
	}

}



bool IncreMarginPerceptron(PointSet &trainPoints, HyperPlane &plane, int dimension, double rho)
{

    //find the max distance between the points and the origin
    double maxDistantce = 0, minDistance = MAX_DOUBLE;
    int n = trainPoints.size();
	//PrintPoints(trainPoints, dimension);
	//compute the distance between two points in d + 1 dimension
	//both points have value 1 in d + 1 dimension
    for (int i = 0; i < n; i++) {
		double curDis = 0;
		for (int j = 0; j < dimension; j++) {
			curDis += trainPoints[i].x[j] * trainPoints[i].x[j];
		}
		curDis = sqrt(curDis + 1);
        if (maxDistantce < curDis) maxDistantce = curDis;
        if (minDistance > curDis) minDistance = curDis;
    }
    cout << "Min: " << minDistance << "  Max: " << maxDistantce << endl;
    double y_guess = maxDistantce;
    while(true) {
		printf("current y_guess:");
        cout << y_guess << endl;
        plane.Clear();
		bool re = BallMarginPerceptron(trainPoints, plane, dimension, y_guess, maxDistantce, minDistance);
        if (re) {
            break;
        } else{
			y_guess = (1 - rho)*y_guess;
            //cout << y_guess << endl;
        }
        if(y_guess < 0.1 ) {////need to consider carefully
            puts("It is non-linear separable!");
            return false;
        }
    }
    return true;
}
