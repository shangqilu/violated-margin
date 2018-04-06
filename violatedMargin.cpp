#include "violatedMargin.h"


//return a hyper-plane violated no more than (1+epsilon)*k points
//and the margin is no less than (1- rho)optimal_margin
//with high successful probability no less than 1 - delta
bool ApproximateViolatedMargin(PointSet &points, HyperPlane &optimal_plane, int dimension,
                               int k, double epsilon, double rho, double delta, int method)
{
    int n = points.size();
    printf("Original Set size %d ...\n", n);

    //similar to 1+epsilon/3
    double temp = epsilon/log(1+epsilon);
    //sampling probability
    double p = 1.0/(temp*log(temp)-temp+1)*(dimension*log(n)-LogFactorial(dimension)+log(1/delta))/k;
    cout << "sample probability p: " << p << endl;

    //make sure in the sample set the plane
    //violating no more than k_prime points
    int k_prime = temp * k * p;
    cout << "k: " << k << " k_prime: " << k_prime << endl;
    if (k_prime > k)
    {
        k_prime = k;
    }
    cout << "k: " << k << " k_prime: " << k_prime << endl;
	PointIndex subSetIndex = Sampling(n, p);
    
	printf("n_prime: %d ...\n", subSetIndex.size());

	bool found = ViolatedMargin(points, subSetIndex, optimal_plane, dimension, k, k_prime, rho, delta, method);

    if (found)
    {
        double margin = 0;
        int real_k = 0;
        bool k_separable = MinimumViolatedDistance(points, optimal_plane, margin, (1+epsilon)*k, real_k);
        if (k_separable)
        {
            puts("Separating original set correctly");
            printf("real_k: %d, margin: %lf, error rate: %lf\n", real_k, margin, 1.0*real_k/n);
            return true;
        }
        else
        {
            return false;
        }
    } else {
        return false;
    }
}




bool ViolatedMargin(PointSet &points, PointIndex &subSetIndex, HyperPlane &optimal_plane, int dimension,
		int original_k, int k, double rho, double delta, int method)
{
    //repeating k^d * exp^d / d^d ln(1/delta) times

    long long times = (long long)(pow(k, dimension) * exp(dimension)
                                /pow(1.0*dimension, dimension) * log(1/delta));
    printf("need repeating times: %lld \n", times);
    double max_margin = 0;

    bool flag = false;
	PointSet coresetDirections;
	PointSet classifyDirections;
	if (method == 2){
		//compute those directions set the angle with epsilon
		//compute direction set for constructing CoreSet
		double alpha = 1.0 / (dimension*(4 * dimension + 1));
		double single_algle = sqrt(rho*alpha / 4);
		double radius = sqrt(dimension) + 1;

		double *angles = new double[dimension - 1];
		puts("computing coreset directions...");
		ComputingDirections(coresetDirections, angles, 1, dimension, 0.2, radius);
		printf("there are %d directions", coresetDirections.size());

		ComputingDirections(classifyDirections, angles, 1, dimension, 0.2, 1);
		printf("there are %d directions", classifyDirections.size());
		delete[]angles;
	}


    for (int curtime = 0; curtime < 5000; curtime++)
    {
        printf("Current repeating time: %d\n", curtime);
        double p = 1.0*dimension/k;
		PointIndex subsubIndex = Sampling(subSetIndex, p);
		PointSet subpoints;
		for (int i = 0; i < subsubIndex.size(); i++) {
			subpoints.push_back(points[subsubIndex[i]]);
		}
		int n = subpoints.size();
        printf("Sub Set size %d ...\n", n);
		if (n < 2) {
			continue;
		}
        HyperPlane plane(dimension);
        bool found = MarginClasification(subpoints, plane, dimension, rho, method, coresetDirections, classifyDirections);
        if (found)
        {
            //test whether separable in the sub set
            double cur_dis = 0;
            int real_k = 0;
			bool k_separable = MinimumSubsetViolatedDistance(points,subSetIndex, plane, cur_dis, k, real_k);
            if (k_separable)
            {
				printf("a solution in the sub set \n");
				bool correct = MinimumViolatedDistance(points, plane, cur_dis, original_k, real_k);
				if (correct && cur_dis > max_margin) {
					flag = true;
					max_margin = cur_dis;
					//optimal_plane = plane;
					printf("a solution in the original set ");
					cout << max_margin << endl;
					PrintHyperPlane(plane);
					CopyHyperPlane(optimal_plane, plane);
				}  
            }
        }
    }
    if (flag)
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

PointSet Sampling(PointSet &points, int dimension, double p)
{
    //srand(time(NULL));   //in reality we should use time seed
    if (p > 1)
    {
        puts("p>1");
        return points;
    }
    int n = points.size();
    PointSet newPoints;
    for (int i = 0; i < n; i++)
    {
        double r = ((double)rand())/RAND_MAX;
        //printf("%d %lf\n", i, r);
        if (r < p)
        {
            newPoints.push_back(points[i]);
        }
    }
    return newPoints;
}


bool MarginClasification(PointSet &points, HyperPlane &plane, int dimension, double rho, 
	int method, PointSet &coresetDirections, PointSet &classifyDirections)
{
    if (method == 0)
    {
        return IncreMarginPerceptron(points, plane, dimension, rho);
    }
    else if (method == 1)
    {
        return LPclassification(points, plane, dimension);
    }
    else if (method == 2)
    {
		return DirectionalWidth(points, plane, dimension, rho, coresetDirections, classifyDirections);
    }
    else
    {
        puts("choose method from 0 to 2");
        return false;
    }
}
