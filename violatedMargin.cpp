#include "violatedMargin.h"
#include "tools.h"
#include "globalVar.h"
#include "perceptron.h"
#include "simplex.h"
#include "directionalWidth.h"

bool ApproximationViolatedMargin(PointSet &points, HyperPlane &optimal_plane)
{
    int n = points.size();
    printf("Original Set size %d ...\n", n);
	
	/*
    *	Similar to k_prime = 1 + epsilon/3 in the book.
	*	when k_prime = (epsilon/log(1+epsilon)) * k * p,
	*	the sampling probability is smaller in theory.
	*/
    double temp = Epsilon/log(1+Epsilon);
    //sampling probability
    double p = 1.0/(temp*log(temp)-temp+1)*(Dim*log(n)-LogFactorial(Dim)+log(1/Delta))/K;
    cout << "sample probability p: " << p << endl;

    //in the sample set the separation plane
    //shoud violate no more than k_prime points
    int k_prime = temp * K * p;
    cout << "k: " << K << " k_prime: " << k_prime << endl;
    if (k_prime > K)
    {
        k_prime = K;
    }
    cout << "k: " << K << " k_prime: " << k_prime << endl;
	PointIndex subSetIndex = Sampling(n, p);
    
	printf("n_prime: %d ...\n", subSetIndex.size());

	bool found = ViolatedMargin(points, subSetIndex, optimal_plane, k_prime);

    if (found)
    {
		return true;
    } else {
        return false;
    }
}




bool ViolatedMargin(PointSet &points, PointIndex &subSetIndex, HyperPlane &optimal_plane, int k_prime)
{
    //should repeat k^d * exp^d / d^d ln(1/delta) times in theory
	//the failing probability delta is set to 0.5 
	long long times = (long long)(pow(k_prime, Dim) * exp(Dim)
					/ pow(1.0*Dim, Dim) * log(1 / Delta));
    printf("need repeating times in theory: %lld \n", times);
    
	double max_margin = 0;
	int real_k;
	int realMaxIterations = 500000;
	int curMax = 10000;
	int find_cnt = 0;
	int max_cnt = 10;
    bool found_solution = false;

	PointSet coresetDirections;
	PointSet classifyDirections;
	if (Method == 2) {
		//compute those directions set the angle with epsilon
		//compute direction set for constructing CoreSet
		double alpha = 1.0 / (Dim*(4 * Dim + 1));
		double single_algle = sqrt(Rho*alpha / 4);
		double radius = sqrt(Dim) + 1;

		double *angles = new double[Dim - 1];
		puts("computing coreset directions...");
		ComputingDirections(coresetDirections, angles, 1, 0.2, radius, 2*M_PI);
		printf("there are %d coreset directions\n", coresetDirections.size());

		ComputingDirections(classifyDirections, angles, 1, 0.2, 1, M_PI);
		printf("there are %d classifying directions\n", classifyDirections.size());
		delete[]angles;
	}





	for (int curtime = 0; curtime < curMax; curtime++)
	{
		if (curtime >= realMaxIterations) {
			break;
		}
		printf("Current repeating time: %d\n", curtime);
		double p = 1.0  / k_prime; //
		if (k_prime * 1.0 / K < 0.8) 
		{
			p = p * Dim;
		}
		PointIndex subsubIndex = Sampling(subSetIndex, p);
		PointSet subsubPoints;
		for (int i = 0; i < subsubIndex.size(); i++) {
			subsubPoints.push_back(points[subsubIndex[i]]);
		}
		int subsubn = subsubPoints.size();
		printf("Sub Set size %d ...\n", subsubn);
		if (subsubn < 2) {
			continue;
		}
        HyperPlane plane(Dim);
		bool found = MarginClasification(subsubPoints, plane, coresetDirections, classifyDirections);
        if (found)
        {
            //test whether separable in the sub set
            double cur_dis = 0;
            real_k = 0;
            int permitted_k = k_prime;
            if (k_prime == K) permitted_k = (1 + Epsilon) * k_prime;
			bool k_separable = MinimumSubsetViolatedDistance(points,subSetIndex, plane, cur_dis, permitted_k, real_k);
            if (k_separable)
            {
				//if we can find a solution in the subset easily,  
				//increase the maxiterations for a better solution.
				curMax = curMax * 2;
				if (curMax > realMaxIterations) curMax = realMaxIterations;
				//printf("find a solution in the sub set \n");
				//test whether separable in the original set
				bool correct = MinimumViolatedDistance(points, plane, cur_dis, (1+Epsilon)*K, real_k);
				if (correct && cur_dis > max_margin) {
					find_cnt++;
					if (find_cnt >= max_cnt) break;
					found_solution = true;
					max_margin = cur_dis;
					//optimal_plane = plane;
					printf("a solution in the original set with margin: ");
					cout << max_margin << endl;
					PrintHyperPlane(plane);
					CopyHyperPlane(optimal_plane, plane);
					break;
				}  
            }
        }
    }

	if (Method == 2) {
		//delete points
		int num = coresetDirections.size();
		for (int i = 0; i < num; i++) {
			delete[]coresetDirections[i]->x;
			delete(coresetDirections[i]);
			coresetDirections[i] = NULL;
		}
		num = classifyDirections.size();
		for (int i = 0; i < num; i++) {
			delete[]classifyDirections[i]->x;
			delete(classifyDirections[i]);
			classifyDirections[i] = NULL;
		}

	}
	if (found_solution)
    {
        return true;
    }
    else
    {
        return false;
    }
}


PointIndex Sampling(int n, double p)
{
	PointIndex index_points;
	if (p > 1) {
		index_points.resize(n);
		for (int i = 0; i < n; i++) {
			index_points[i] = i;
		}
		return index_points;
	}
	for (int i = 0; i < n; i++)
	{
		double r = ((double)rand()) / RAND_MAX;
		//printf("%d %lf\n", i, r);
		if (r < p)
		{
			index_points.push_back(i);
		}
	}
	return index_points;

}

PointIndex Sampling(PointIndex &points, double p)
{
	PointIndex index_points;
	if (p > 1) {
		return points;
	}
	int n = points.size();
	for (int i = 0; i < n; i++)
	{
		double r = ((double)rand()) / RAND_MAX;
		//printf("%d %lf\n", i, r);
		if (r < p)
		{
			index_points.push_back(points[i]);
		}
	}
	return index_points;

}




inline bool MarginClasification(PointSet &points, HyperPlane &plane,
				PointSet &coresetDirections, PointSet &classifyDirections)
{
    if (Method == 0)
    {
		return LPclassification(points, plane);
    }
    else if (Method == 1)
    {
		return IncreMarginPerceptron(points, plane);
    }
	else if (Method == 2) {
		return DirectionalWidth(points, plane, coresetDirections, classifyDirections);
	}
	else
    {
        puts("should choose method from 0 to 2");
        return false;
    }
}
