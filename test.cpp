#include "test.h"

using namespace std;


void TestSimplex()
{
    /*
    char train_data[] = "data/separable_test_2/titanic_train_data.asc";
    char train_label[] = "data/separable_test_2/titanic_train_label.asc";
    */
    char filename[] = "data/iris_4.txt";
    int dimension = 4;
    PointSet train_points = LoadDataLibSVMFormat(filename, dimension);
    cout << train_points.size() << endl;
    //test for Simplex
    HyperPlane plane = HyperPlane(dimension);
    bool found = LPclassification(train_points, plane, dimension);
    if (found) {
        puts("find a solution");
    }
    int n = train_points.size();
    double cur_distance = 0;
    bool separable = MinimumSeparableDistance(train_points, plane, cur_distance);
    if (separable) {
        PrintHyperPlane(plane, plane.d);
        cout << "cur_dis: " << cur_distance << endl;
    }

}

void TestSimplex2()
{
    double *input_A[3];
    double A[3][3]= {{1,1,3},{2,2,5},{4,1,2}};

    for(int i = 0; i < 3; i++)
        input_A[i] = A[i];
    double input_B[3] = {30,24,36};
    double input_C[3] = {3,1,2};


    /*

    double *input_A[2];
    double A[2][2]= {{2, -1},{1,-5}};
    double input_B[2] = {2,-4};
    double input_C[2] = {2,-1};
    */

    LPresult result = Simplex(input_A, input_B, input_C, 3, 3);
    PrintLPresult(result, 3);
}

void TestPerceptron()
{
    /*
    char train_data[] = "data/separable_test_2/titanic_train_data.asc";
    char train_label[] = "data/separable_test_2/titanic_train_label.asc";
    int dimension = 2;
    PointSet train_points = LoadData(train_data, train_label, dimension);
    */
    //test for Margin Perceptron
    //SimplePerceptron(train_points, dimension);

    char filename[] = "data/iris_4.txt";
    int dimension = 4;
    PointSet train_points = LoadDataLibSVMFormat(filename, dimension);

    HyperPlane plane = HyperPlane(dimension);
    bool found = IncreMarginPerceptron(train_points, plane, dimension, 0.3);
    if(found) {
        puts("find a solution");
        double cur_distance = 0;
        bool separable = MinimumSeparableDistance(train_points, plane, cur_distance);
        if (separable) {
            PrintHyperPlane(plane, plane.d);
            cout << "cur_dis: " << cur_distance << endl;
        }
    }else {
        puts("there is no solution");
    }

}

void TestGaussianEquation()
{
    //test for guassian_equiation
    double a[3][3] = { 1,1,1,0,4,-1,2,-2,1 };
    double **A = new double*[3];
    for (int i = 0; i < 3; i++)
    {
        A[i] = a[i];
    }
    double b[3] = { 6,5,1};
    double x[3];
    GaussianEquation(A,b,x,3);
    for (int i = 0; i < 3; i++ )
        cout << x[i] << " ";
    cout << endl;

}

void TestGaussianInverse()
{
    //test for gaussian inverse
    double a[3][3] = { 1,0,1,1,1,2,3,4,2};
    double **A = new double*[3], **B = new double*[3];
    for (int i = 0; i < 3; i++)
    {
        A[i] = a[i];
        B[i] = new double[3];
    }
    GaussianInverseMatrix(A,B,3);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            cout << B[i][j] << " ";
        }
        cout << endl;
    }
}

void TestBoudingBox()
{
    char train_data[] = "data/separable_test_2/titanic_train_data.asc";
    char train_label[] = "data/separable_test_2/titanic_train_label.asc";
    int dimension = 2;
    PointSet train_points = LoadData(train_data, train_label, dimension);

    MinimumBoudingBox(train_points, dimension);

}


void TestSimpleCoreSet()
{
    char train_data[] = "data/separable_test_2/titanic_train_data.asc";
    char train_label[] = "data/separable_test_2/titanic_train_label.asc";
    int dimension = 2;
    PointSet points = LoadData(train_data, train_label, dimension);
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
    PointSet simpleset = SimpleCoreSet(points, MainTainT, dimension, 0.1);
    PrintPoints(simpleset, dimension);
}

void TestSmallerCoreSetandComputingDirection()
{
    char train_data[] = "data/separable_test_2/titanic_train_data.asc";
    char train_label[] = "data/separable_test_2/titanic_train_label.asc";
    int dimension = 2;
    PointSet points = LoadData(train_data, train_label, dimension);

    double epsilon = 0.1;
    double alpha = 1.0/(dimension*(4*dimension+1));
    double delta = sqrt(epsilon*alpha/4);
    double radius = sqrt(dimension) + 1;


    double *angles = new double[dimension-1];
    PointSet directionPoints;
    ComputingDirections(directionPoints, angles, 1, dimension, delta, radius);
    printf("there are %d points in all directions\n", directionPoints.size());
    //PrintPoints(directionPoints, dimension);

    PointSet smallerSet = SmallerCoreSet(points, directionPoints, dimension, epsilon);
    puts("after transform");
    PrintPoints(smallerSet, dimension);
}

void TestComputingDirections()
{
    PointSet directionPoints;
    int dimension = 4;
    double *angles = new double[dimension];
    double epsilon = 0.1;
    ComputingDirections(directionPoints, angles, 1, dimension, epsilon, 1);
    printf("there are %d directions\n", directionPoints.size());
}

void TestOneDimensionClassification()
{
    char train_data[] = "data/separable_test_2/titanic_train_data.asc";
    char train_label[] = "data/separable_test_2/titanic_train_label.asc";
    int dimension = 2;
    PointSet points = LoadData(train_data, train_label, dimension);

    double epsilon = 0.001;
    double radius = 1;

    double *angles = new double[dimension-1];
    PointSet directionPoints;
    ComputingDirections(directionPoints, angles, 1, dimension, epsilon, radius);
    printf("there are %d points in all directions\n", directionPoints.size());

    int n = points.size();
    HyperPlane optimal_plane = HyperPlane(dimension);
    double max_margin = 0;
    //
    for (int i = 0; i < directionPoints.size(); i++)
    {
        HyperPlane cur_plane = HyperPlane(dimension);

        bool separable = OneDimensionClassification(points, cur_plane, directionPoints[i], radius);
        //compute the margin of current plane
        if (separable)
        {
            //check each hyperplane direction with the original data
            double cur_distance = MAX_DOUBLE;
            for (int i = 0; i < n; i++)
            {
                double tmp_dis = Distance(cur_plane, points[i], dimension);
                if (tmp_dis < cur_distance)
                {
                    cur_distance = tmp_dis;
                }
            }
            if (cur_distance > max_margin)
            {
                PrintHyperPlane(cur_plane, dimension);
                cout << "cur_dis: " << cur_distance << endl;
                optimal_plane = cur_plane;
                max_margin = cur_distance;
            }
        }
    }
}

void TestDirectionalWidth()
{
    /*
    char train_data[] = "data/separable_test_2/titanic_train_data.asc";
    char train_label[] = "data/separable_test_2/titanic_train_label.asc";
    int dimension = 2;
    char train_data[] = "data/titanic_3/titanic_train_data.asc";
    char train_label[] = "data/titanic_3/titanic_train_label.asc";
    int dimension = 3;
    */
    //PointSet points = LoadData(train_data, train_label, dimension);

    double epsilon = 0.1;

    char filename[] = "data/skin_nonskin_3.txt";
    int dimension = 3;
    PointSet train_points = LoadDataLibSVMFormat(filename, dimension);

    HyperPlane plane = HyperPlane(dimension);
    bool found = DirectionalWidth(train_points, plane, dimension, epsilon);
    if(found) {
        puts("find a solution");
        double cur_distance = 0;
        bool separable = MinimumSeparableDistance(train_points, plane, cur_distance);
        if (separable) {
            PrintHyperPlane(plane, plane.d);
            cout << "cur_dis: " << cur_distance << endl;
        }

    }else {
        puts("there is no solution");
    }

}


void TestViolatedMargin(int method)
{
    char filename[] = "data/skin_nonskin_3.txt";
    //char filename[] = "data/svm_guide_4.txt";
    int dimension = 3;
    PointSet trainPoints = LoadDataLibSVMFormat(filename, dimension);
    HyperPlane plane = HyperPlane(dimension);
    int k = trainPoints.size() * 0.06;
    double epsilon = 0.1;
    double rho = 0.1;
    double delta = 0.5;
    ApproximateViolatedMargin(trainPoints, plane, dimension, k, epsilon, rho, delta, method);
}
