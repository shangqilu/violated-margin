#include "headers.h"
#include "tools.h"
#include "globalVar.h"
#include "violatedMargin.h"
#include "dataGenerator.h"
#include <errno.h>

using namespace std;

static char *line = NULL;
static int max_line_len;

static char* readline(FILE *input)
{
	if (fgets(line, max_line_len, input) == NULL)
		return NULL;

	while (strrchr(line, '\n') == NULL)
	{
		max_line_len *= 2;
		line = (char *)realloc(line, max_line_len);
		int len = (int)strlen(line);
		if (fgets(line + len, max_line_len - len, input) == NULL)
			break;
	}
	return line;
}

void exit_input_error(int line_num)
{
	printf("Wrong input format at line %d\n", line_num);
	exit(1);
}
void LoadDataLibSVMFormat(int n, int dimension, char* filename, PointSet &points)
{
	FILE *fp = fopen(filename, "r");
	if (fp == NULL) {
		fprintf(stderr, "can't open input file %s\n", filename);
		exit(1);
	}
	max_line_len = 1024;
	line = (char *)malloc(max_line_len * sizeof(char));
	int n_instance = 0;
	int n_features = 0;
	char *c_label, *endptr;
	char *idx, *val;
	bool positive_label = false;
	bool negative_label = false;
	clock_t startreadtime = clock();

	while (readline(fp) != NULL){
		n_instance++;
		if (n_instance % 10000 == 0) {
			printf("%d\n", n_instance);
		}
		int inst_max_index = 0;
		c_label = strtok(line, " \t\n");
		if (c_label == NULL)
		{
			exit_input_error(n_instance);
		}

		errno = 0;
		int cur_label = (int)strtol(c_label, &endptr, 10);
		if (endptr == c_label || errno != 0 || *endptr != '\0')
			exit_input_error(n_instance);

		if (cur_label != -1 && cur_label != 1) {
			printf("label should be -1 or 1 !\n");
			exit_input_error(n_instance);
		}
		//cout << cur_label << endl;
		double *x = new double[dimension];
		Point *curPt = new Point(dimension, x);
		if (n_instance > 245057) {
			cout << n_instance << endl;
		}
		if (cur_label == 1) {
			positive_label = true;
			curPt->y = 1;
		}
		else {
			negative_label = true;
			curPt->y = -1;
		}
		
		while (1)
		{
			idx = strtok(NULL, ":");
			val = strtok(NULL, " \t");

			if (val == NULL)
			{
				break;
			}
				

			errno = 0;
			int cur_index = (int)strtol(idx, &endptr, 10);
			if (endptr == idx || errno != 0 || *endptr != '\0' || cur_index <= inst_max_index
				|| cur_index > dimension)
			{
				printf("index must start from 1 to d in ascending order\n");
				exit_input_error(n_instance);
			}	
			else
			{
				inst_max_index = cur_index;
			}

			errno = 0;
			double cur_val = strtod(val, &endptr);
			if (endptr == val || errno != 0 || (*endptr != '\0' && !isspace(*endptr)))
			{
				exit_input_error(n_instance);
			}
			curPt->x[cur_index - 1] = cur_val;
			//cout << cur_index << " " << cur_val << " " << n_instance << endl;
		}
		points.push_back(curPt);
		if (n_features < inst_max_index) {
			n_features = inst_max_index;
		}
	}

	if (n_instance != n) {
		PrintError("The number of instances in the file is not equal to n");
	} 
	if (positive_label == false || negative_label == false) {
		PrintError("There is only one class of label.");
	}

	clock_t endreadtime = clock();
	double duration = (endreadtime - startreadtime) / (double)CLOCKS_PER_SEC;
	printf("read data time: %lf\n", duration);
	
	free(line);
	fclose(fp);
}

void Usage()
{
	printf(
	"k violation margin algorithm\n"
	"Options:\n"
	"-n {integer} the number of points\n"
	"-d {integer} the dimensionlity of points\n"
	"-p {double}  the largest classification error rate accepted\n" 
	"			  or the percent of noise points\n"
	"-f {string}  the full path of training file name\n"
	"			  or generated dataset file name\n"
	"-e {double}  epsilon: the approximation for error rate (default 0.1)\n"
	"-r {double}  rho: the approximation for margin (default 0.1)\n"
	"-m {integer}  the choosen method (default 0)\n"
	"   0:  the linear programming algorithm\n"
	"   1:  the perceptron algorithm\n" 
	"   2:  the directional width algorithm\n"
	"   3:  generating data points\n"
	);

}

void WritePlanetoFile(char *filename, HyperPlane &plane, double margin, int real_k, double error_rate)
{
	FILE *fp = fopen(filename, "w");
	if (fp == NULL) {
		printf("cannot open file %s\n", filename);
		exit(1);
	}
	fprintf(fp, "Dim:%d\n", plane.d);
	for (int i = 0; i < plane.d; i++) {
		fprintf(fp, "%lf ", plane.w[i]);
	}
	fprintf(fp, "%lf\n", plane.b);
	fprintf(fp, "margin: %lf, reak_k: %d, error_rate: %lf\n", margin, real_k, error_rate);
	fclose(fp);
}


void InitializeGlobalVariables(int n, int d, double p, int k, double epsilon, double rho, int method)
{
	N = n;
	Dim = d;
	Percent = p;
	K = k;
	Epsilon = epsilon;
	Rho = rho;
	Method = method;
}



int main(int argc, char **argv)
{
#ifndef __DEBUG__
	srand(time(NULL));
#endif // __DEBUG__

	srand(time(NULL));
	//for train
	

	int n = -1;
	int d = -1;
	double p = -1;
	double epsilon = 0.1;
	double rho = 0.1;
	int method = 0;
	char file_name[200]; 
	char model_file_name[200];
	bool flag = true;

	//parse command line
	int i;
	for (i = 1; i < argc; i++) {
		if (argv[i][0] != '-') {
			flag = false;
			break;
		}
		if (++i > argc) {
			flag = false;
			break;
		}
		switch(argv[i-1][1]) {
			case 'n' :{
				n = atoi(argv[i]);
				if (n <= 0) flag = false;
				printf("n: %d\n", n);
				break;
			}
			case 'd' :{
				d = atoi(argv[i]);
				if (d <= 0) flag = false;
				printf("d: %d\n", d);
				break;
			}
			case 'p' :{
				p = atof(argv[i]);
				if (p < 0 || p > 1) flag = false;
				printf("p: %lf\n", p);
				break;
			}
			case 'f' :{
				strcpy(file_name, argv[i]);
				char *model = strrchr(argv[i], '/');
				if (model == NULL) {
					model = argv[i];
				}else {
					++ model;
				}
				sprintf(model_file_name, "%s.model", model);
				break;
			}
			case 'e' :{
				epsilon = atof(argv[i]);
				if (epsilon < 0) flag = false;
				break;
			}
			case 'r' :{
				rho = atof(argv[i]);
				if (rho < 0 || rho > 1) flag = false;
				break;
			}
			case 'm' :{
				method = atoi(argv[i]);
				if (method < 0 || method > 3) flag = false;
				printf("m: %d\n", method);
				break;
			}
			default: {
				flag = false;
				break;
			}

		}

	}
	if (n <= 0 || d <= 0 || p < 0)
	{
		flag = false;
		printf("=========================================\n");
		printf("You should indicate -n, -d, -p, -f at least.\n");
	}
	if (!flag) {
		printf("=========================================\n");
		printf("There is something wrong in your input.\n");
		printf("=========================================\n");
		Usage();
		exit(1);
	}
	


	//given n, d, k, train_datafilename, epsilon, rho
	int k = (int) n * p;
	if (k < 0) {
		k = 0;
	}
	else if (k > n) {
		k = n;
	}

	//initialize global variables
	InitializeGlobalVariables(n, d, p, k, epsilon, rho, method);
	
	if (method == 3) {
		printf("Generating points\n");
		GenMarginDataSet(file_name, 1.0 / d, 1, n, k);
		return 0;
	}


	//start time
	clock_t start = clock();

	//create PointSet
	PointSet points;
	points.reserve(n);

	//read data
	LoadDataLibSVMFormat(n, d, file_name, points);
	//PrintPoints(points, d);
	cout << points.size() << endl;

	//call the approximation margin algorithm
	HyperPlane plane(Dim);
	bool found = ApproximationViolatedMargin(points, plane);
	
	if (found) {
		double margin;
		int real_k;
		MinimumViolatedDistance(points, plane, margin, (1 + Epsilon) * K, real_k);
		PrintHyperPlane(plane);
		printf("margin: %lf, reak_k: %d, error_rate: %lf\n", margin, real_k, real_k*1.0 / N);
		WritePlanetoFile(model_file_name, plane, margin, real_k, real_k*1.0 / N);
	}
	else {
		printf("Failed to find a solution");
	}


	puts("Deleting points...");
	//release PointSet
	for (int i = 0; i < N; i++) {
		//if (i > 0 && i % 10000 == 0) printf("%d\n", i);
		delete[] points[i]->x;
		delete (points[i]);
		points[i] = NULL;
	}
	points.clear();
	//end time
	clock_t endtime = clock();


	//write result into a file/the console 
	double duration = (endtime - start) / (double) CLOCKS_PER_SEC;
	printf("\nprogram ended...\n");
	printf("Traing time: %lfs\n", duration);


	return 0;
}
