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
        Point newpt(dimension);
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
	f_label.close();
    return trainPointSet;
}


//load data from single file
PointSet LoadDataLibSVMFormat(char* filename, int dimension)
{
	PointSet points;
	points.reserve(5000000);
	FILE *fp = fopen(filename, "r");
	if (fp == NULL) {
		PrintError("cannot open file");
	}
	int cnt = 0;
	while (!feof(fp)) {
		if (cnt % 10000 == 0) cout << cnt << endl;
		char c;
		while (1) {
			c = fgetc(fp);
			if (c == EOF) {
				return points;
			}
			if (c == '-' || isdigit(c)) {
				ungetc(c, fp);
				break;
			}
		}
		Point cur_pt(dimension);
		int label;
		fscanf(fp, "%d", &label);
		if (label == 1) {
			cur_pt.y = 1;
		}
		else {
			cur_pt.y = -1;
		}

		int i;
		double content;
		while ((c = fgetc(fp)) != '\n') {
			if (c == EOF) {
				break;
			}
			if (isdigit(c)) {
				ungetc(c, fp);
				fscanf(fp, "%d:%lf", &i, &content);
				//cout << i << " " << content << " " << cnt << endl;
				cur_pt.x[i - 1] = content;
			}
		}
		points.push_back(cur_pt);
		cnt++;
	}
	fclose(fp);
	if (points.size() < points.capacity())
	{
		points.shrink_to_fit();
	}
	return points;
}


PointSet CopyPoints(PointSet &points, int dimension)
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

void CopyHyperPlane(HyperPlane &plane1, HyperPlane &plane2)
{
    for (int i = 0; i < plane1.d; i++)
    {
        plane1.w[i] = plane2.w[i];
    }
    plane1.b = plane2.b;
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

double Dot(Point &pt1, Point &pt2, int dimension)
{
    double sum = 0;
	for (int i = 0; i < dimension; i++)
    {
        sum += pt1.x[i]*pt2.x[i];
    }
    return sum;
}

double Dot(vector<double> x, vector<double> y, int dimension)
{
	double sum = 0;
	for (int i = 0; i < dimension; i++)
	{
		sum += x[i] * y[i];
	}
	return sum;
}

Point PointMinus(Point &pt1, Point &pt2, int dimension)
{  
    Point ans(dimension);
	for (int i = 0; i < dimension; i++)
    {
        ans.x[i] = pt1.x[i] - pt2.x[i];
    }
    return ans;
}

void PrintHyperPlane(HyperPlane &plane)
{
    cout << "w: (";
    for (int i = 0; i < plane.d; i++)
    {
        cout << plane.w[i];
        if (i != plane.d - 1)
        {
            cout << ", ";
        }
    }
    cout << ") ";
    cout << "b: " << plane.b << endl;
}

double Distance(HyperPlane &plane, Point &pt, int dimension)
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

double Distance(Point &pt1, Point &pt2, int dimension)
{
    double dis = 0;
	for (int i = 0; i < dimension; i++)
    {
        dis += (pt1.x[i] - pt2.x[i]) * (pt1.x[i] - pt2.x[i]);
    }
    return sqrt(dis);
}


bool MinimumSeparableDistance(PointSet &points, HyperPlane &plane, double &min_dis)
{
    int dimension = plane.d;
    int n = points.size();
    min_dis = 0;
    //cout << "***" << n << endl;
    double min_margin = MAX_DOUBLE;
    for (int i = 0; i < n; i++)
    {
		if (points[i].y*(Dot(plane.w, points[i].x, dimension) + plane.b) < ZERO_ERROR)
        {
            return false;
        }
        double cur_dis = Distance(plane, points[i], plane.d);
        if (cur_dis < min_margin ) {
            min_margin = cur_dis;
        }
    }
	if (fabs(min_margin - MAX_DOUBLE) < ZERO_ERROR) {
        return false;
    } else {
        min_dis = min_margin;
        return true;
    }
}


bool MinimumViolatedDistance(PointSet &points, HyperPlane &plane, double &min_dis, int k, int &real_k)
{
    int dimension = plane.d;
    int n = points.size();
    min_dis = 0;
    //cout << "***" << n << endl;
    double min_margin = MAX_DOUBLE;
    int wrong_num = 0;
    for (int i = 0; i < n; i++)
    {
		if (points[i].y*(Dot(plane.w, points[i].x, dimension) + plane.b) < ZERO_ERROR)
        {
            wrong_num ++;
        }else {
            double cur_dis = Distance(plane, points[i], plane.d);
            if (cur_dis < min_margin ) {
                min_margin = cur_dis;
            }
        }
    }
    real_k = wrong_num;
    if (wrong_num > k) {
        printf("wrong numbers : %d\n", wrong_num);
        return false;
    }
	if (fabs(min_margin - MAX_DOUBLE) < ZERO_ERROR) {
        return false;
    } else {
        min_dis = min_margin;
        return true;
    }
}


bool MinimumSubsetViolatedDistance(PointSet &points, PointIndex &index, HyperPlane &plane, double &min_dis, int k, int &real_k)
{
	int dimension = plane.d;
	int n = index.size();
	min_dis = 0;
	//cout << "***" << n << endl;
	double min_margin = MAX_DOUBLE;
	int wrong_num = 0;
	for (int i = 0; i < n; i++)
	{
		if (points[index[i]].y*(Dot(plane.w, points[index[i]].x, dimension) + plane.b) < ZERO_ERROR)
		{
			wrong_num++;
		}
		else {
			double cur_dis = Distance(plane, points[index[i]], plane.d);
			if (cur_dis < min_margin) {
				min_margin = cur_dis;
			}
		}
	}
	real_k = wrong_num;
	if (wrong_num > k) {
		printf("wrong numbers : %d\n", wrong_num);
		return false;
	}
	if (fabs(min_margin - MAX_DOUBLE) < ZERO_ERROR) {
		return false;
	}
	else {
		min_dis = min_margin;
		return true;
	}
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
		if (fabs(colMax) < ZERO_ERROR)
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
        delete []M[i];
    }
    delete []M;
    return true;
}

bool GaussianInverseMatrix(double** A, double**B, int n)
{
    double **Matrix = new double*[n];
    double **Inverse = new double*[n];
    //Copy A and Initial B unit matrix
    for (int i = 0; i < n; i++)
    {
        Matrix[i] = new double[n];
        Inverse[i] = new double[n];
        for (int j = 0; j < n; j++)
        {
            Matrix[i][j] = A[i][j];
        }
        for (int j = 0; j < n; j++)
        {
            Inverse[i][j] = (i==j)? 1:0;
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
		if (fabs(colMax) < ZERO_ERROR)
        {
            //puts("There is no solution to the inverse matrix");
            //singular matrix
            for (int i = 0; i < n; i++)
            {
                delete []Matrix[i];
                delete []Inverse[i];
            }
            delete []Matrix;
            delete []Inverse;
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

                temp = Inverse[i][j];
                Inverse[i][j] = Inverse[maxLineIndex][j];
                Inverse[maxLineIndex][j] = temp;
            }
        }

        // set Matrix[i][i] = 1
        double cur = Matrix[i][i];
        for (int j = 0; j < n; j++)
        {
            Matrix[i][j] /= cur;
            //cout <<  Matrix[i][j] << " "<<  endl;
            Inverse[i][j] /= cur;
            //cout << Inverse[i][j] << " " <<endl;
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
                    Inverse[j][k] -= Inverse[i][k] * temp;
                }
            }
        }

    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            B[i][j] = Inverse[i][j];
        }
    }

    //puts("*****2");
    //delete memory
    for (int i = 0; i < n; i++)
    {
        delete []Matrix[i];
        delete []Inverse[i];
    }
    delete []Matrix;
    delete []Inverse;
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


void PrintPoints(PointSet &points, int dimension)
{
    for (int i = 0; i < points.size(); i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            cout << points[i].x[j] << " ";
        }
        cout << "label: " << points[i].y;
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
    delete []x;
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
        delete []T[i];
    }
    delete []T;
}


int LargestPrimeBelow(int m)
{
    if (m < 2) {
        printf("there is no prime under %d", m);
        return 2;
    }
    int cur = m;
    while(cur --)
    {
        if (is_prime(cur))
        {
            //printf("%d's largest prime is %d\n", m, cur);
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


long long Factorial(int d)
{
    long long ans = 1;
    for (int i = 2; i <= d; i++)
    {
        ans = ans * i;
    }
    return ans;
}

double LogFactorial(int d)
{
    double ans = 0;
    for (int i = 1; i <= d; i++)
    {
        ans += log(i);
    }
    return ans;
}


void PrintError(string msg)
{
	printf("%s\n", msg.c_str());
	exit(0);
}