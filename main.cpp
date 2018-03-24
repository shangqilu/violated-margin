#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<cstring>
#include<algorithm>
#include<sstream>
#include<fstream>
#include<string>
#include<vector>
using namespace std;


struct HyperPlane
{
    double *w;
    double b;
    HyperPlane(int d)
    {
        this->w = new double[d];
        memset(this->w, 0, sizeof(this->w));
        this->b = 0;
    }
};

struct Point
{
    double *x;
    int y;
    Point(int d)
    {
        this->x = new double[d];
    }
    Point(int d, double *x, int y)
    {
        this->x = new double[d];
        memcpy(this->x, x, sizeof(x));
        this->y = y;
    }
};

typedef vector<Point> PointSet;
PointSet train_points;

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
        while(line >> tmp) {
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
    while(f_label >> cur_label) {
        trainPointSet[index].y = cur_label;
        index++;
    }
    printf("Loading %d training examples.\n", index);

    return trainPointSet;
}


double dot(double* w, double *x, int dimension)
{
    double sum = 0;
    for (int i = 0; i < dimension; i++) {
        sum += (*w) * (*x);
        w++;
        x++;
    }
    return sum;
}

void Perceptron(PointSet trainPoints, int dimension)
{
    int n = trainPoints.size();
    HyperPlane plane = HyperPlane(dimension);
    int iter_cnt = 0;
    while(true)
    {
        printf("Iteration %d\n", ++iter_cnt);
        int i;
        for(i = 0; i < n; i++)
        {
            Point cur_pt = trainPoints[i];
            if(cur_pt.y * (dot(plane.w, cur_pt.x, dimension) + plane.b)   <= 0)
            {
                for (int j = 0; j < dimension; j++)
                {
                    plane.w[j] += cur_pt.y*cur_pt.x[j];
                }
                plane.b += cur_pt.y;
                printf(" %d  ", i);
                break;
            }
        }
        if(i == n)
            break;
    }
}


int main()
{
    char train_data[] = "data/titanic_3/titanic_train_data.asc";
    char train_label[] = "data/titanic_3/titanic_train_label.asc";
    int dimension = 3;
    train_points = LoadData(train_data, train_label, dimension);

    Perceptron(train_points, dimension);



    return 0;
}
