#include "perceptron.h"
#include "headers.h"
#include "simplex.h"
#include <iostream>

using namespace std;


PointSet LoadData(char* filename, char* label_filename, int dimension)
{
    PointSet trainPointSet;
    ifstream f_data, f_label;
    f_data.open(filename);
    string input;

    //loading the training data
    while(getline(f_data, input))
    {
        stringstream line(input);
        Point newpt = Point(dimension);
        double tmp = 0;
        int index = 0;
        while(line >> tmp)
        {
            newpt.x[index] = tmp;
            index ++;
        }
        trainPointSet.push_back(newpt);
    }
    f_data.close();
    //loading the training label
    f_label.open(label_filename);
    int cur_label;
    int index = 0;
    while(f_label >> cur_label)
    {
        trainPointSet[index].y = cur_label;
        index++;
    }
    printf("Loading %d training examples.\n", index);

    return trainPointSet;
}






int main()
{


    char train_data[] = "data/separable_test_2/titanic_train_data.asc";
    char train_label[] = "data/separable_test_2/titanic_train_label.asc";
    int dimension = 2;
    PointSet train_points = LoadData(train_data, train_label, dimension);

    //test for Simplex
    HyperPlane plane = LPclassification(train_points, dimension, 2);
    PrintHyperPlane(plane, plane.d);

    //test for Margin Perceptron
    //SimplePerceptron(train_points, dimension);

    //IncreMarginPerceptrion(train_points, dimension, 0.3);

    //char lpfile[] = "data/LPtest/simplelp.asc";

    /*
    double *input_A[3];
    double A[3][3]= {{1,1,3},{2,2,5},{4,1,2}};

    for(int i = 0; i < 3; i++)
        input_A[i] = A[i];
    double input_B[3] = {30,24,36};
    double input_C[3] = {3,1,2};



    for(int i = 0; i < 3; i++)
        input_A[i] = A[i];

    double *input_A[2];
    double A[2][2]= {{2, -1},{1,-5}};
    double input_B[2] = {2,-4};
    double input_C[2] = {2,-1};


    LPresult result = Simplex(input_A, input_B, input_C, 3, 3);
    PrintLPresult(result, 3);

    */

    return 0;
}
