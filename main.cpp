#include "test.h"
#include "headers.h"
#include <iostream>

using namespace std;


//sampling with probability p
//copy those points into a new Point Set
PointSet Sampling(PointSet points, int dimension, double p)
{
    //srand(time(NULL));   //in reality we should use time seed
    if (p > 1) {
        puts("p>1");
        return points;
    }
    int n = points.size();
    PointSet newPoints;
    for (int i = 0; i < n; i++)
    {
        double r = ((double)rand())/RAND_MAX;
        if (r < p)
        {
            newPoints.push_back(points[i]);
        }
    }
    return newPoints;
}


bool MarginClasification(PointSet points, HyperPlane &plane, int dimension, double rho, int method)
{
    if (method == 0) {
        return IncreMarginPerceptron(points, plane, dimension, rho);
    }else if (method == 1)
    {
        return LPclassification(points, plane, dimension);
    }else if (method == 2){
        return DirectionalWidth(points, plane, dimension, rho);
    } else {
        puts("choose method from 0 to 2");
        return false;
    }
}

//return a hyper-plane violated no more than (1+epsilon)*k points
//and the margin is no less than (1- rho)optimal_margin
//with high successful probability no less than 1 - delta
bool ViolatedMargin(PointSet points, HyperPlane &optimal_plane, int dimension, int k, double epsilon, double rho, double delta, int method)
{
    int n = points.size();
    printf("Original Set size %d ...\n", n);
    //similar to 1+epsilon/3
    double temp = epsilon/log(1+epsilon);
    //the first time sampling probability
    double p = 1.0/(temp*log(temp)-temp+1)*(dimension*log(n)-LogFactorial(dimension)+log(1/delta))/k;
    cout << "p: " << p << endl;
    //make sure in the sample set the plane
    //violating no more than k_prime points
    double k_prime = temp * k * p;
    cout << "k_prime: " << k_prime << endl;
    //srand(time(NULL));
    PointSet subpoints = Sampling(points, dimension, p);
    printf("Sub Set size %d ...\n", subpoints.size());

    //the second time sampling probability : dimension/k_prime
    //repeating k_prime^d * exp^d / d^d ln(1/delta) times

    long long int times = (long long)(pow(k_prime, dimension) * exp(dimension)/pow(1.0*dimension, dimension) * log(1/delta));
    cout << times << endl;
    printf("need repeating times: %lld \n", times);
    double max_margin = 0;
    for (int curtime = 0; curtime < 1000000; curtime++)
    {
        printf("Current repeating time: %d\n", curtime);
        double p = 1.0*dimension/k_prime;
        PointSet subsubpoints = Sampling(subpoints, dimension, p);
        printf("Sub sub Set size %d ...\n", subsubpoints.size());
        HyperPlane plane = HyperPlane(dimension);
        bool found = MarginClasification(subsubpoints, plane, dimension, rho, method);
        if (found) {
            puts("find a solution in the subsub set");
            double cur_dis = MinimumSeparableDistance(subpoints, plane);
            if (cur_dis > max_margin) {
                max_margin = cur_dis;
                //optimal_plane = plane;
                printf("a solution on the sub set ");
                cout << max_margin << endl;
                PrintHyperPlane(plane, dimension);
                CopyHyperPlane(optimal_plane, plane);
                break;
            }
        }
    }

    if (max_margin > 0) {
        PrintHyperPlane(optimal_plane, dimension);
        puts("find a hyperplane classifying the subpoints correctly");
    }

    double real_margin = MinimumSeparableDistance(points, optimal_plane);
    if (real_margin > 0) {
        puts("Separating original set correctly");
        return true;
    }else {
        return false;
    }
}




int main()
{


    //TestSimplex();

    //TestOneDimensionClassification();
    //TestPerceptron();

    //TestSimpleCoreSet();

    char filename[] = "data/skin_nonskin_3.txt";
    int dimension = 3;
    PointSet trainPoints = LoadDataLibSVMFormat(filename, dimension);
    HyperPlane plane = HyperPlane(dimension);
    int k = trainPoints.size() * 0.01;
    double epsilon = 0.1;
    double rho = 0.1;
    double delta = 0.5;
    double method = 1;
    ViolatedMargin(trainPoints, plane, dimension, k, epsilon, rho, delta, method);

    //TestDirectionalWidth();



    return 0;
}
