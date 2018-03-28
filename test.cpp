#include "test.h"

using namespace std;


void TestSimplex()
{
    char train_data[] = "data/separable_test_2/titanic_train_data.asc";
    char train_label[] = "data/separable_test_2/titanic_train_label.asc";
    int dimension = 2;
    PointSet train_points = LoadData(train_data, train_label, dimension);

    //test for Simplex
    HyperPlane plane = LPclassification(train_points, dimension, 2);
    PrintHyperPlane(plane, plane.d);
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
    char train_data[] = "data/separable_test_2/titanic_train_data.asc";
    char train_label[] = "data/separable_test_2/titanic_train_label.asc";
    int dimension = 2;
    PointSet train_points = LoadData(train_data, train_label, dimension);

    //test for Margin Perceptron
    //SimplePerceptron(train_points, dimension);

    HyperPlane plane = IncreMarginPerceptron(train_points, dimension, 0.3);
    PrintHyperPlane(plane, plane.d);
}

void TestGaussianEquation()
{
    //test for guassian_equiation
    double a[3][3] = { 1,1,1,0,4,-1,2,-2,1 };
    double **A = new double*[3];
    for (int i = 0; i < 3; i++) {
        A[i] = a[i];
    }
    double b[3] = { 6,5,1};
    double x[3];
    GaussianEquation(A,b,x,3);
    for (int i = 0; i < 3; i++ ) cout << x[i] << " ";
    cout << endl;

}

void TestGaussianInverse()
{
    //test for gaussian inverse
    double a[3][3] = { 1,0,1,1,1,2,3,4,2};
    double **A = new double*[3], **B = new double*[3];
    for (int i = 0; i < 3; i++) {
        A[i] = a[i];
        B[i] = new double[3];
    }
    GaussianInverseMatrix(A,B,3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
        {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    for (int i = 0; i < 3; i++) {
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
