#include "perceptron.h"
#include "headers.h"
#include "simplex.h"
#include "test.h"
#include <iostream>

using namespace std;




int main()
{


    char train_data[] = "data/separable_test_2/titanic_train_data.asc";
    char train_label[] = "data/separable_test_2/titanic_train_label.asc";
    int dimension = 2;
    PointSet train_points = LoadData(train_data, train_label, dimension);
    TestBoudingBox();



    return 0;
}
