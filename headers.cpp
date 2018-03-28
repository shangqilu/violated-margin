#include "headers.h"

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

PointSet CopyPoints(PointSet points, int dimension)
{
    int n = points.size();
    PointSet newPoints;
    for (int i = 0; i < n; i++) {
        Point tmp = Point(dimension);
        for (int j = 0; j < dimension; j++)
        {
            tmp.x[i] = points[i].x[i];
        }
        tmp.y = points[i].y;
        newPoints.push_back(tmp);
    }
    return newPoints;
}

double Dot(double* w, double *x, int dimension)
{
    double sum = 0;
    for (int i = 0; i < dimension; i++)
    {
        sum += (*w) * (*x);
        w++;
        x++;
    }
    return sum;
}

double Dot(Point pt1, Point pt2, int dimension)
{
    double sum = 0;
    for (int i = 0; i < dimension; i++)
    {
        sum += pt1.x[i]*pt2.x[i];
    }
    return sum;
}

Point PointMinus(Point pt1, Point pt2, int dimension)
{
    Point ans = Point(dimension);
    for (int i = 0; i < dimension; i++)
    {
        ans.x[i] = pt1.x[i] - pt2.x[i];
    }
    return ans;
}

void PrintHyperPlane(HyperPlane plane, int dimension)
{
    cout << "w: (";
    for (int i = 0; i < dimension; i++)
    {
        cout << plane.w[i];
        if (i != dimension - 1)
        {
            cout << ", ";
        }
    }
    cout << ") ";
    cout << "b: " << plane.b << endl;
}

double Distance(HyperPlane plane, Point pt, int dimension)
{
    double dis = 0, denominator = 0;
    for (int i = 0; i < dimension; i++)
    {
        dis += plane.w[i]*pt.x[i];
        denominator += plane.w[i]*plane.w[i];
    }
    dis += plane.b;
    return fabs(dis)/sqrt(denominator);
}

double Distance(Point pt1, Point pt2, int dimension)
{
    double dis = 0;
    for(int i = 0; i < dimension; i++)
    {
        dis += (pt1.x[i] - pt2.x[i]) * (pt1.x[i] - pt2.x[i]);
    }
    return sqrt(dis);
}



bool GaussianEquation(double** A, double* b, double* x, int n)
{


    double** M = new double*[n];
    for (int i = 0; i < n; i++)
    {
        M[i] = new double[n+1];
        for (int j = 0; j < n; j++)
        {
            M[i][j] = A[i][j];
        }
        M[i][n] = b[i];
    }
    for (int i = 0; i < n; i++)
    {
        //find the pivot
        double colMax = fabs(M[i][i]);
        int maxLineIndex = i;
        for(int j = i + 1; j < n; j++)
        {
            //find the largest element in current column
            if(fabs(M[j][i]) > colMax)
            {
                colMax = fabs(M[j][i]);
                maxLineIndex = j;
            }
        }
        if(fabs(colMax) < ZERO)
        {
            puts("There is no solution to the equation set");
            //singular matrix
            for (int i = 0; i < n; i++)
            {
                delete M[i];
            }
            delete M;
            return false;
        }
        if (maxLineIndex != i)
        {
            //swap line k and maxLineIndex
            for (int j = 0; j < n + 1; j++)
            {
                double temp = M[i][j];
                M[i][j] = M[maxLineIndex][j];
                M[maxLineIndex][j] = temp;
            }
        }


        //eliminate cur column below row k
        for(int j = i + 1; j < n; j++)
        {
            for (int k = i + 1; k < n + 1; k++)
            {
                M[j][k] = M[j][k] - M[i][k] * M[j][i] / M[i][i];
            }
        }
    }
    //calculate the solution
    for (int i = n - 1; i >= 0; i--)
    {
        x[i] = M[i][n];
        for (int j = i + 1; j < n; j++)
        {
            x[i] -= M[i][j]*x[j];
        }
        x[i] /= M[i][i];
    }
    for (int i=0; i<n; i++)
    {
        delete M[i];
    }
    delete M;
    return true;
}

bool GaussianInverseMatrix(double** A, double**B, int n)
{
    double **Matrix = new double*[n];
    //Copy A and Initial B unit matrix
    for (int i = 0; i < n; i++)
    {
        Matrix[i] = new double[n];
        for (int j = 0; j < n; j++)
        {
            Matrix[i][j] = A[i][j];
        }
        for (int j = 0; j < n; j++)
        {
            B[i][j] = (i==j)? 1:0;
        }
    }
    for (int i = 0; i < n; i++)
    {
        //find the pivot
        double colMax = fabs(Matrix[i][i]);
        int maxLineIndex = i;
        for(int j = i + 1; j < n; j++)
        {
            //find the largest element in current column
            if(fabs(Matrix[j][i]) > colMax)
            {
                colMax = fabs(Matrix[j][i]);
                maxLineIndex = j;
            }
        }
        if(fabs(colMax) < ZERO)
        {
            puts("There is no solution to the inverse matrix");
            //singular matrix
            for (int i = 0; i < n; i++)
            {
                delete Matrix[i];
            }
            delete Matrix;
            return false;
        }
        if (maxLineIndex != i)
        {
            //swap line k and maxLineIndex
            for (int j = 0; j < n; j++)
            {
                double temp = Matrix[i][j];
                Matrix[i][j] = Matrix[maxLineIndex][j];
                Matrix[maxLineIndex][j] = temp;

                temp = B[i][j];
                B[i][j] = B[maxLineIndex][j];
                B[maxLineIndex][j] = temp;
            }
        }

        // set Matrix[i][i] = 1
        double cur = Matrix[i][i];
        for (int j = 0; j < n; j++)
        {
            Matrix[i][j] /= cur;
            //cout <<  Matrix[i][j] << " "<<  endl;
            B[i][j] /= cur;
            //cout << B[i][j] << " " <<endl;
        }
        //puts("*****2");
        for (int j = 0; j < n; j++)
        {
            if (j != i)
            {
                double temp = Matrix[j][i];
                for (int k = 0; k < n; k++)
                {
                    Matrix[j][k] -= Matrix[i][k] * temp;
                    B[j][k] -= B[i][k] * temp;
                }
            }
        }

    }
    //puts("*****2");
    //delete memory
    for (int i = 0; i < n; i++)
    {
        delete Matrix[i];
    }
    delete Matrix;

    return true;
}




void PrintMatrix(double **T, int dimension)
{
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            cout << T[i][j] << " ";
        }
        cout << endl;
    }
}

void TransformingPoints(PointSet &points, double **T, int point_dimension)
{
    int n = points.size();
    double *x = new double[point_dimension+1];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < point_dimension + 1; j++)
        {
            x[j] = 0;
            for (int k = 0; k < point_dimension + 1; k++)
            {
                if (k < point_dimension)
                    x[j] += points[i].x[k] * T[j][k];
                else
                    x[j] += T[j][k];
            }
        }
        for (int j = 0; j < point_dimension; j++)
        {
            points[i].x[j] = x[j];
        }
    }
    delete x;
}

void MatrixMultiply(double **A, double **B, double **C, int dimension)
{
    double **T = new double*[dimension];
    for (int i = 0; i < dimension; i++)
    {
        T[i] = new double[dimension];
    }

    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            T[i][j] = 0;
            for (int k = 0; k < dimension; k++)
            {
                T[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            C[i][j] = T[i][j];
        }
    }
    for (int i = 0; i < dimension; i++)
    {
        delete T[i];
    }
    delete T;
}


int LargestPrimeBelow(int m)
{
    if (m < 2) {
        puts("there is no prime under m");
        return 0;
    }
    int cur = m;
    while(cur --)
    {
        if (is_prime(cur))
        {
            printf("%d's largest prime is %d\n", m, cur);
            return cur;
        }
    }
}

bool is_prime(int n)
{
    long long i;
    for (i = 2; i * i <= n; i++) {
        if ((n%i) == 0) return false;
    }
    return true;
}

