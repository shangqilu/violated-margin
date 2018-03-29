#include "test.h"
#include "headers.h"
#include <iostream>

using namespace std;


//sampling with probability p
//copy those points into a new Point Set
PointSet Sampling(PointSet points, int dimension, double p)
{
    //srand(time(NULL));   //in reality we should use time seed
    int n = points.size();
    PointSet newPoints;
    for (int i = 0; i < n; i++)
    {
        double r = ((double)rand())/RAND_MAX;
        if (r < p)
        {
            Point tmp = Point(dimension);
            for (int j = 0; j < dimension; j++)
            {
                tmp.x[j] = points[i].x[j];
            }
            tmp.y = points[i].y;
            newPoints.push_back(tmp);
        }
    }
    return newPoints;
}

//return a hyper-plane violated no more than (1+epsilon)*k points
//and the margin is no less than (1- rho)optimal_margin
//with high successful probability no less than 1 - delta
bool ViolatedMargin(PointSet points, HyperPlane &plane, int dimension, int k, double epsilon, double rho, double delta)
{
    int n = points.size();
    //similar to 1+epsilon/3
    double temp = epsilon/log(1+epsilon);
    //the first time sampling probability
    double p = 1.0/(temp*log(temp)-temp+1)*(dimension*log(n)-LogFactorial(dimension)+log(1/delta))/k;

    //make sure in the sample set the plane
    //violating no more than k_prime points
    double k_prime = temp * k * p;
    //srand(time(NULL));
    PointSet newPoints;
    if (p > 1)
    {
        newPoints = points;
    }
    else     //do sampling
    {
        newPoints = Sampling(points, dimension, p);
    }


    //the second time sampling probability : dimension/k_prime
    //repeating k_prime^d * exp^d / d^d ln(1/delta) times
    int times = ceil(pow(k_prime, dimension) * exp(dimension)/pow(1.0*dimension, dimension) * log(1/delta));
    for (int curtime = 0; curtime < times; curtime++)
    {
        double p = dimension/k_prime;
        PointSet subPoints = Sampling(newPoints, dimension, p);
        //bool found = LPclassification();
        //bool found = IncreMarginPerceptron(,);
        //bool found = DirectionalWidth(subPoints,dimension, epsilon);

    }


}




int main()
{


    //TestOneDimensionClassification();

    //TestOneDimensionClassification();
    //TestPerceptron();

    TestDirectionalWidth();

    return 0;
}
