#include "perceptron.h"
#include "tools.h"
#include "globalVar.h"
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
		for (i = 0; i < n; i++)
		{
			Point *cur_pt = trainPoints[i];
			if (cur_pt->y * (Dot(plane.w, cur_pt->x, dimension) + plane.b) <= 0)
			{
				for (int j = 0; j < dimension; j++)
				{
					plane.w[j] += cur_pt->y*cur_pt->x[j];
				}
				plane.b += cur_pt->y;
				//printf("break %d\n", i);
				break;
			}
		}
		if (i == n)
			break;
	}
	if (iter_cnt == MaxIterations) {
		puts("It is non-linearly separable!");
	}
	else PrintHyperPlane(plane);
}


bool BallMarginPerceptron(PointList &pointlist, HyperPlane &plane, double y_guess)
{
	/*
	*	compute the min distance from a d dimension circle in d + 1 dimension to origin
	*	in the last dimension it has value 1
	*/
	double R = 0;
	PointList::iterator iter;
	for (iter = pointlist.begin(); iter != pointlist.end(); ++iter) {
		double d_dis = 0;
		for (int j = 0; j < Dim; j++)
		{
			d_dis += (*iter)->x[j] * (*iter)->x[j];
		}
		d_dis = sqrt(d_dis);
		double tmp_max = sqrt((d_dis + y_guess) * (d_dis + y_guess) + 1);
		if (tmp_max > R){
			R = tmp_max;
		}
	}
	//R^2/y_min^2 cannot be the upper bound
	//since there are infinite number of points
	double y_min = R / (2 * (Dim + 1));
	//1/dimension is a large margin
	//we assume data has a large margin
	long long maxIters = (int)((R*R) / (y_min*y_min*y_guess) + 0.5);
	int iter_cnt = 0;	
	//maxIters *= 100; //its a parameter
	if (maxIters > MaxIterations) {
		maxIters = MaxIterations;
	}
	//printf("maxIters: %d\n", maxIters);
	//cout << "R: " << R << " y_guess: " << y_guess << " y_min: " << y_min << endl;

	for (iter_cnt = 0; iter_cnt < maxIters; iter_cnt++)
	{
		//printf("Iteration %d\n", iter_cnt);
		iter = pointlist.begin();
		bool update = false;
		for (; iter != pointlist.end(); ++iter) 
		{
			if ((*iter)->y * (Dot(plane.w, (*iter)->x, Dim) + plane.b) <= 0)
			{
				for (int j = 0; j < Dim; j++)
				{
					plane.w[j] += (*iter)->y*(*iter)->x[j];
				}
				plane.b += (*iter)->y;
				update = true;
			}
			//check whether (d+1) dimension cur plane passes thorough 
			//d dimension Ball(points[i], y_guess);
			else if (Distance(plane, *(*iter), Dim) < y_guess) {
				double denominator = 0;
				for (int j = 0; j < Dim; j++) {
					denominator += plane.w[j] * plane.w[j];
				}
				if (denominator < ZERO_ERROR) {
					puts("wrong");
					continue;
				}
				//W_t+1 = W_t+ y*(x - y*y_guess*W_t/||W_t||)
				for (int j = 0; j < Dim; j++)
				{
					plane.w[j] += (*iter)->y*((*iter)->x[j] -
						(*iter)->y * y_guess * plane.w[j] / sqrt(denominator));
				}
				plane.b += (*iter)->y;
				update = true;
			}
			if (update) {
				pointlist.push_front(*iter);
				pointlist.erase(iter);
				break;
			}
		}
		if (!update)
		{
			break;
		}
	}
	if (iter_cnt >= maxIters) {
		//puts("It is non-linearly separable based on y_guess");
		return false;
	}
	else {
		//PrintHyperPlane(plane);
		return true;
	}

}
/*
*	compute the min distance from a d dimension circle in d + 1 dimension
*	in the last dimension it has value 1
*/

bool IncreMarginPerceptron(PointSet &trainPoints, HyperPlane &plane)
{
	PointList pointlist;
	for (int i = 0; i < trainPoints.size(); i++) {
		pointlist.push_back(trainPoints[i]);
	}
    //find the max distance between the points and the origin
    double maxDistantce = 0;
    int n = trainPoints.size();

	/*
	*	compute the largest Distance from a point to the a origin
	*	and guess the margin starting from it. 
	*/
    for (int i = 0; i < n; i++) {
		double curDis = 0;
		for (int j = 0; j < Dim; j++) {
			curDis += trainPoints[i]->x[j] * trainPoints[i]->x[j];
		}
        if (maxDistantce < sqrt(curDis)) maxDistantce = sqrt(curDis);
    }
    //cout << "Max: " << maxDistantce << endl;
    double y_guess = maxDistantce;
    while(true) {

        plane.Clear();
		bool re = BallMarginPerceptron(pointlist, plane, y_guess);
        if (re) {
            break;
        } else{
			y_guess = (1 - Rho)*y_guess;
        }
        if(y_guess < 1.0/(2*Dim) ) {////need to consider carefully
            //puts("data set is non-linear separable!");
            return false;
        }
    }
    return true;
}
