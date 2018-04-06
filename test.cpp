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
    //PrintPoints(train_points, dimension);
    //test for Simplex
    HyperPlane plane(dimension);
    bool found = LPclassification(train_points, plane, dimension);
	puts("begin check");
    if (found) {
        puts("find a solution");
    } else {
        puts("no solution");
    }
    int n = train_points.size();
    double cur_distance = 0;
    bool separable = MinimumSeparableDistance(train_points, plane, cur_distance);
    if (separable) {
        PrintHyperPlane(plane);
        cout << "cur_dis: " << cur_distance << endl;
    }

}

void TestSimplex2()
{
    /*

    double *input_A[3];
    double A[3][3]= {{1,1,3},{2,2,5},{4,1,2}};

    for(int i = 0; i < 3; i++)
        input_A[i] = A[i];
    double input_B[3] = {30,24,36};
    double input_C[3] = {3,1,2};
    */    

    double *input_A[2];
    double A[2][2]= {{2, -1},{1,-5}};
    for(int i = 0; i < 2; i++)
        input_A[i] = A[i];
    double input_B[2] = {2,-4};
    double input_C[2] = {2,-1};
    
    LPresult result = Simplex(input_A, input_B, input_C, 2, 2);
    if (result.flag == 1) {
        puts("find a solution");
        PrintLPresult(result, 2);
        puts("***");
    } else {
        puts("no solution");
    }
    //
}



void TestPerceptron()
{
    
    /*char train_data[] = "data/separable_test_2/titanic_train_data.asc";
    char train_label[] = "data/separable_test_2/titanic_train_label.asc";
    int dimension = 2;
    PointSet train_points = LoadData(train_data, train_label, dimension);
	PrintPoints(train_points, dimension);*/

    //test for Margin Perceptron
    //SimplePerceptron(train_points, dimension);
    
    char filename[] = "data/iris_4.txt";
    int dimension = 4;
	PointSet train_points = LoadDataLibSVMFormat(filename, dimension);
    //SimplePerceptron(train_points, dimension);
	
    HyperPlane plane(dimension);
    //bool found = MarginPerceptron(train_points, plane, dimension, 0.1, 3, 0.3);
    
    bool found = IncreMarginPerceptron(train_points, plane, dimension, 0.1);
    if(found) {
        puts("find a solution");
        double cur_distance = 0;
		PrintHyperPlane(plane);
        bool separable = MinimumSeparableDistance(train_points, plane, cur_distance);
        if (separable) {
            PrintHyperPlane(plane);
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
    for (int i = 0; i < 3; i++) {
        delete []B[i];
    }
    delete []A;
    delete []B;
}


void TestBoundingBox()
{
    char train_data[] = "data/separable_test_2/titanic_train_data.asc";
    char train_label[] = "data/separable_test_2/titanic_train_label.asc";
    int dimension = 2;
    PointSet points = LoadData(train_data, train_label, dimension);
	PrintPoints(points, dimension);
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

    BoundingBox box(dimension);
    RecursionMinimumBoudingBox(points, box, MainTainT, dimension, dimension);
    //PrintBoudingBox(box);

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
    for (int i = 0; i < dimension+1; i++)
    {
        delete []MainTainT[i];
    }
    delete []MainTainT;

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
	puts("****");
    PrintPoints(simpleset, dimension);
    for (int i = 0; i < dimension+1; i++)
    {
        delete []MainTainT[i];
    }
    delete []MainTainT;
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
	//169 directions
    ComputingDirections(directionPoints, angles, 1, dimension, delta, radius);
    printf("there are %d points in all directions\n", directionPoints.size());
    //PrintPoints(directionPoints, dimension);

    PointSet smallerSet = SmallerCoreSet(points, directionPoints, dimension, epsilon);
    puts("after transform");
	//(0,0) is deleted
    PrintPoints(smallerSet, dimension);
    delete []angles;
}

void TestComputingDirections()
{
    PointSet directionPoints;
    int dimension = 3;

	double epsilon = 0.1;
	double alpha = 1.0 / (dimension*(4 * dimension + 1));
	double delta = sqrt(epsilon*alpha / 4);
	double radius = sqrt(dimension) + 1;

    double *angles = new double[dimension];
    //double epsilon = 0.1;
	ComputingDirections(directionPoints, angles, 1, dimension, delta, radius);
	printf("there are %d points in all directions\n", directionPoints.size());
    delete []angles;
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
    HyperPlane optimal_plane(dimension);
    double max_margin = 0;
    //
    for (int i = 0; i < directionPoints.size(); i++)
    {
        HyperPlane cur_plane(dimension);

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
                PrintHyperPlane(cur_plane);
                cout << "cur_dis: " << cur_distance << endl;
				CopyHyperPlane(optimal_plane, cur_plane);
                max_margin = cur_distance;
            }
        }
    }
    delete []angles;
}

void TestDirectionalWidth()
{
    
    char train_data[] = "data/separable_test_2/titanic_train_data.asc";
    char train_label[] = "data/separable_test_2/titanic_train_label.asc";
    int dimension = 2;
	double epsilon = 0.1;
	PointSet points = LoadData(train_data, train_label, dimension);
	/*
    //
    char filename[] = "data/iris_4.txt";
    int dimension = 4;
    PointSet train_points = LoadDataLibSVMFormat(filename, dimension);
	*/
    HyperPlane plane(dimension);
	bool found = DirectionalWidth(points, plane, dimension, epsilon);
    if(found) {
        puts("find a solution");
        double cur_distance = 0;
		bool separable = MinimumSeparableDistance(points, plane, cur_distance);
        if (separable) {
            PrintHyperPlane(plane);
            cout << "cur_dis: " << cur_distance << endl;
        }

    }else {
        puts("there is no solution");
    }

}


void TestViolatedMargin(int method)
{
    char filename[] = "D:/My work/ConsoleApplication1/ConsoleApplication1/data/Margin_datasetE6D5R5.txt";
	//char filename[] = "data/poker.t";
    int dimension = 5;
	PointSet trainPoints = LoadDataLibSVMFormat(filename, dimension);
    HyperPlane plane(dimension);
    int k = trainPoints.size() * 0.01;
    double epsilon = 0.1;
    double rho = 0.1;
    double delta = 0.5;
    //PrintPoints(trainPoints, dimension);
    puts("******");
    ApproximateViolatedMargin(trainPoints, plane, dimension, k, epsilon, rho, delta, method);
}

void TestSampling()
{
    char filename[] = "D:/code/violated-margin/data/skin_nonskin_3.txt";
    int dimension = 3;
    PointSet trainPoints = LoadDataLibSVMFormat(filename, dimension);
    //PrintPoints(trainPoints, dimension);
    PointSet subset = Sampling(trainPoints, dimension, 0.5);
    //puts("*****");
    //PrintPoints(subset, dimension);
    cout << subset.size();
}



void TestDataGenerator()
{
	char filename[] = "data/Margin_datasetE6D5R5.txt";

	int totalNum = 1000000;
	int noiseNum = totalNum * 0.01;
	int dimension = 5;
	double margin = 1.0 / dimension;
	double radius = 5;
	GenMarginDataSet(filename, dimension, margin, radius, totalNum, noiseNum);

}