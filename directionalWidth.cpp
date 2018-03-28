#include "directionalWidth.h"
#include <ctime>


using namespace std;

void TwoApproxiDiameter(PointSet points, int dimension, int &s, int &t)
{
    //compute the 2-approximation diameter
    int n = points.size();
    //srand(time(NULL));
    while(1)
    {
        s = rand()% n;
        double max_dis = -MAX_DOUBLE;
        t = -1;
        for (int i = 0; i < n; i++)
        {
            if (i != s)
            {
                double cur_dis = Distance(points[s], points[i], dimension);
                if (cur_dis > max_dis)
                {
                    t = i;
                    max_dis = cur_dis;
                }
            }
        }
        int flag = 0;
        Point direction = PointMinus(points[s], points[t], dimension);
        Point axis = Point(dimension);
        axis.x[dimension-1] = 1;
        //the hyper plane s and t defined is not parallel to last dimension coordinate axis
        if (fabs(Dot(direction, axis, dimension)) > ZERO)
        {
            break;
        }
    }
}




void RecursionMinimumBoudingBox(PointSet points, BoudingBox &curbox, double **MainTainT, int dimension, int realDimension)
{
    if (dimension == 0)
    {
        return;
    }
    puts("Computing the approximate minimum bouding box recursively");
    printf("cur dimension %d...\n", dimension);

    int n = points.size();
    if (n < 2)
    {
        puts("there are no more than 2 points!");
        return;
    }
    int s, t;
    TwoApproxiDiameter(points, dimension, s, t);
    printf("current s and t: %d %d\n", s, t);
    //all points is above point t
    Point d_axis = PointMinus(points[s], points[t], dimension);
    //find those d-1 basis

    double **CurT = new double*[realDimension+1];
    //store those d basis vectors
    for (int i = 0; i < realDimension+1; i++)
    {
        CurT[i] = new double[realDimension+1];
        for (int j = 0; j < realDimension+1; j++)
        {
            CurT[i][j] = 0;
        }
    }
    //filling the transforming matrix
    for (int i = dimension; i < realDimension+1; i++)
    {
        CurT[i][i] = 1;
    }

    for (int i = 0; i < dimension; i++)
    {
        CurT[i][dimension-1] = d_axis.x[i];
        CurT[i][realDimension] = points[t].x[i];
    }

    //find those d-1 basis
    double **A, *b, *x;
    if (dimension > 1)
    {
        A = new double*[dimension-1];
        for (int i = 0; i < dimension-1; i++)
        {
            A[i] = new double[dimension-1];
        }
        b = new double[dimension-1];
        x = new double[dimension-1];
    }
    for (int number_variables = 1; number_variables <= dimension-1; number_variables++)
    {
        // there are number_variables variables
        //set the previous dimension-i variables to 1
        //form the equation set
        for (int i = 0; i < number_variables; i++)
        {
            double sum = 0;
            for (int j = 0; j < dimension - number_variables; j++)
            {
                sum += CurT[j][dimension-1-i]*(1-points[t].x[j]);
            }
            for (int j = dimension - number_variables; j < dimension; j++)
            {
                A[i][j- (dimension-number_variables)] = CurT[j][dimension-1-i];
                sum -= CurT[j][dimension-1-i] * points[t].x[j];
            }
            b[i] = -sum;
        }
        bool solve = GaussianEquation(A, b, x, number_variables);
        if (!solve)
        {
            puts("the equation set is insouble!");
            break;
        }

        for (int j = 0; j < dimension - number_variables; j++)
        {
            CurT[j][dimension-1-number_variables] = 1 - points[t].x[j];
        }
        for (int j = dimension - number_variables; j < dimension; j++)
        {
            CurT[j][dimension-1-number_variables] = x[j-(dimension-number_variables)] - points[t].x[j];
        }
    }


    //Print Current basis
    puts("current basis ...");
    PrintMatrix(CurT, realDimension+1);

    //compute the inverse Matrix of cur_T
    double **Inverse_M = new double*[realDimension+1];
    //store those d basis vectors
    for (int i = 0; i < realDimension+1; i++)
    {
        Inverse_M[i] = new double[realDimension+1];
    }
    bool re = GaussianInverseMatrix(CurT, Inverse_M, realDimension+1);
    if (!re)
    {
        puts("there is no inverse Matrix about CurT");
    }
    puts("inverse Matrix...");
    PrintMatrix(Inverse_M, realDimension+1);


    TransformingPoints(points, Inverse_M, realDimension); ////!!!!!!

    MatrixMultiply(Inverse_M, MainTainT, MainTainT, realDimension + 1);
    puts("Maintain Matrix...");
    PrintMatrix(MainTainT, realDimension+1);

    //find the point with the largest value in current dimension
    double maxValue = 0;
    double max_index = t;
    for (int i = 0; i < n; i++)
    {
        if (points[i].x[dimension-1] > maxValue)
        {
            maxValue = points[i].x[dimension-1];
            max_index = i;
        }
    }
    curbox.L[dimension-1] = points[t].x[dimension-1];
    curbox.U[dimension-1] = points[max_index].x[dimension-1];
    for (int i = 0; i < realDimension+1; i++)
    {
        delete CurT[i];
        delete Inverse_M[i];
    }
    delete CurT;
    delete Inverse_M;

    RecursionMinimumBoudingBox(points, curbox, MainTainT, dimension-1, realDimension);


}

void MinimumBoudingBox(PointSet points, int dimension)
{
    double **MainTainT = new double*[dimension+1];
    //store those d basis vectors
    for (int i = 0; i < dimension+1; i++)
    {
        MainTainT[i] = new double[dimension+1];
        for (int j = 0; j < dimension+1; j++)
        {
            MainTainT[i][j] = 0;
        }
    }
    for (int i = 0; i < dimension + 1; i++)
        MainTainT[i][i] = 1;
    BoudingBox box = BoudingBox(dimension);
    RecursionMinimumBoudingBox(points, box, MainTainT, dimension, dimension);
    PrintBoudingBox(box);

    puts("the new points");
    int n = points.size();
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j  < dimension; j++)
        {
            cout << points[i].x[j] << " ";
        }
        cout << endl;
    }
    PrintMatrix(MainTainT, dimension+1);
}


void PrintBoudingBox(BoudingBox box)
{
    for (int i = 0; i < box.d; i++)
    {
        printf("%d dimension: %lf %lf\n", i, box.L[i], box.U[i]);
    }
}

PointSet SimpleCoreSet(PointSet points, int dimension, double epsilon)
{

    BoudingBox box = BoudingBox(dimension);
    double **MainTainT = new double*[dimension+1];
    //store those d basis vectors
    for (int i = 0; i < dimension+1; i++)
    {
        MainTainT[i] = new double[dimension+1];
        for (int j = 0; j < dimension+1; j++)
        {
            MainTainT[i][j] = 0;
        }
    }
    for (int i = 0; i < dimension + 1; i++)
    {
        MainTainT[i][i] = 1;
    }

    RecursionMinimumBoudingBox(points, box, MainTainT, dimension, dimension);
    PrintBoudingBox(box);



    //map the bounding box into a hypercube with length 2, centering origin
    double **HypercubeM = new double*[dimension+1];
    for (int i = 0; i < dimension; i++)
    {
        HypercubeM[i] = new double[dimension+1];

        for (int j = 0; j < dimension+1; j++) {
            if(j < dimension) {
                if (j == i) {
                    HypercubeM[i][j] = 2/box.U[j];
                }else{
                    HypercubeM[i][j] = 0;
                }
            }else {
                HypercubeM[i][j] = -1;
            }

        }
    }
    HypercubeM[dimension] = new double[dimension+1];
    for (int j = 0; j < dimension + 1; j++)
    {
        if (j < dimension) {
            HypercubeM[dimension][j] = 0;
            box.U[j] = 1;
            box.L[j] = -1;
        }else{
            HypercubeM[dimension][j] = 1;
        }
    }

    TransformingPoints(points, HypercubeM, dimension);
    MatrixMultiply(HypercubeM, MainTainT, MainTainT, dimension+1);

    double C_d = 1.0/(dimension*(4*dimension+1)); //a parameter
    int M = ceil(4/(epsilon*C_d)); //number of intervals in each dimension
    printf("M: %d \n");
    int n = points.size();
    int prime = LargestPrimeBelow(1.5*n);
    HashTable table = HashTable(1.5*n, prime);

    for (int i = 0; i < n; i++)
    {
        int key = 0;
        for (int j = 0; j < dimension - 1; j++) {
            int cur_k = floor((points[i].x[j] + 1)/(2/M));
            key = (key * M + cur_k) % prime;
        }
        table.Insert(points[i], dimension, i, key);
    }
    vector<int> index_coreset = table.Travel();


}

PointSet SmallerCoreSet(PointSet points)
{

}


HyperPlane DirectionalWidth(PointSet points, double epsilon)
{
    //find +1 points and -1 points

    //compute those directions set the angle with epsilon



    //compute the epsilon-core set respectively according to those directions



    //compute a hyperplane along each direction

    //check each hyperplane direction with the original data


    //return the largest margin hyperplane

}


void ComputingDirections(PointSet &points, int dimension)
{
    if (dimension == 1) {


        return;
    }
}

