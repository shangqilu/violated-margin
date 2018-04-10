#include "test.h"
#include "globalVar.h"
#include "tools.h"
using namespace std;


void TestSimplex()
{
    /*
    char train_data[] = "data/separable_test_2/titanic_train_data.asc";
    char train_label[] = "data/separable_test_2/titanic_train_label.asc";
    */

	

}

void TestSimplex2()
{
    

    double *input_A[3];
    double A[3][3]= {{1,1,3},{2,2,5},{4,1,2}};

    for(int i = 0; i < 3; i++)
        input_A[i] = A[i];
    double input_B[3] = {30,24,36};
    double input_C[3] = {3,1,2};
        

    //double *input_A[2];
    //double A[2][2]= {{2, -1},{1,-5}};
    //for(int i = 0; i < 2; i++)
    //    input_A[i] = A[i];
    //double input_B[2] = {2,-4};
    //double input_C[2] = {2,-1};
    
	/*double *input_A[3];
	double A[3][4] = { { 1, 1, 1, 1 }, { 0.5, -5.5, -2.5, 9 }, { 0.5,-1.5,-0.5,1} };
	double input_B[3] = {1,0,0};
	double input_C[4] = { -1, 7, 1, 2};
	for (int i = 0; i < 3; i++)input_A[i] = A[i];*/

    LPresult result = Simplex(input_A, input_B, input_C, 3, 3);
    if (result.flag == 1) {
        puts("find a solution");
        PrintLPresult(result, 4);
        puts("***");
    } else {
        puts("no solution");
    }
    //
}

void TestSimplexfromFile()
{
	int m = 66;
	int n = 23;
	char filename[] = "data/cycleLP-cycle.txt";
	
	double **A = new double*[m];
	double *B = new double[m];
	double *C = new double[n];
	for (int i = 0; i < m; i++)
	{
		A[i] = new double[n];
	}

	ifstream f_data, f_label;
	f_data.open(filename);
	string input;

	//loading the training data
	int cnt = 0;
	while (getline(f_data, input))
	{
		stringstream line(input);
		double tmp = 0;
		int index = 0;
		if (cnt < m) {
			while (line >> tmp)
			{
				A[cnt][index] = tmp;
				index++;
			}
		}
		else if (cnt == m){
			while (line >> tmp)
			{
				B[index] = tmp;
				index++;
			}
		}
		else {
			while (line >> tmp)
			{
				C[index] = tmp;
				index++;
			}
		}
		cnt++;
	}
	f_data.close();

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
		{
			cout << A[i][j] << " ";
		}
		cout << endl;
	}
	for (int i = 0; i < m; i++) cout << B[i] << " ";
	cout << endl;
	for (int i = 0; i < n; i++) cout << C[i] << ' ';
	cout << endl;
	LPresult result = Simplex(A, B, C, m, n);
	if (result.flag == 1) {
		puts("find a solution");
		PrintLPresult(result, n);
		puts("***");
	}
	else {
		puts("no solution");
	}

	for (int i = 0; i < m; i++)
	{
		delete A[i];
	}
	delete B;
	delete C;
}

