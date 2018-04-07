#include "headers.h"
#include "test.h"
#include <iostream>
#include <errno.h>
using namespace std;



static char *line = NULL;
static int max_line_len;

static char* readline(FILE *input)
{
	int len;

	if (fgets(line, max_line_len, input) == NULL)
		return NULL;

	while (strrchr(line, '\n') == NULL)
	{
		max_line_len *= 2;
		line = (char *)realloc(line, max_line_len);
		len = (int)strlen(line);
		if (fgets(line + len, max_line_len - len, input) == NULL)
			break;
	}
	return line;
}

void exit_input_error(int line_num)
{
	fprintf(stderr, "Wrong input format at line %d\n", line_num);
	exit(1);
}
PointSet LoadData_line(char* filename)
{
	PointSet points;

	//first look through all lines check the data
	//obtain numbers of instance and features

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

	while (readline(fp) != NULL){
		n_instance++;
		int inst_max_index = 0;
		c_label = strtok(line, " \t\n");
		if (c_label == NULL)
		{
			exit_input_error(n_instance);
		}

		errno = 0;
		int tmp = (int)strtol(c_label, &endptr, 10);
		if (endptr == c_label || errno != 0 || *endptr != '\0')
			exit_input_error(n_instance);


		/*if (tmp != 1 && tmp != -1)
		{
		puts("label should be 1 or -1");
		exit_input_error(n_instance);
		}*/

		while (1)
		{
			idx = strtok(NULL, ":");
			val = strtok(NULL, " \t");

			if (val == NULL)
				break;

			errno = 0;
			int cur_index = (int)strtol(idx, &endptr, 10);
			if (endptr == idx || errno != 0 || *endptr != '\0' || cur_index <= inst_max_index)
				exit_input_error(n_instance);
			else
				inst_max_index = cur_index;

			errno = 0;
			double cur_val = strtod(val, &endptr);
			if (endptr == val || errno != 0 || (*endptr != '\0' && !isspace(*endptr)))
				exit_input_error(n_instance);
		}
		if (n_features < inst_max_index) {
			n_features = inst_max_index;
		}
	}
	rewind(fp);
	cout << n_instance << " ";
	cout << n_features << endl;
	points.reserve(n_instance);

	for (int i = 0; i < n_instance; i++)
	{
		//Point cur_pt(n_features);
		readline(fp);
		c_label = strtok(line, " \t\n");

		int label = (int)strtol(c_label, &endptr, 10);

		/*if (label == 1) {
		cur_pt.y = 1;
		}
		else {
		cur_pt.y = -1;
		}*/
		while (1)
		{
			idx = strtok(NULL, ":");
			val = strtok(NULL, " \t");

			if (val == NULL)
				break;

			int cur_index = (int)strtol(idx, &endptr, 10);
			double cur_val = strtod(val, &endptr);

			cout << cur_index << " " << cur_val << " " << i << endl;
			//cur_pt.x[cur_index-1] = cur_val;
		}
		//points.push_back(cur_pt);
	}

	free(line);
	fclose(fp);
	return points;
}




int main(int nargs, char **args)
{
#ifndef __DEBUG__
	srand(time(NULL));
#endif // __DEBUG__
	srand(time(NULL));
	//TestDirectionalWidth();
	int method = 0;
	//TestViolatedMargin(method);
	///TestSimplexfromFile();
	//TestSimplex();
	//TestPerceptron();
	//TestDataGenerator();
	//TestDirectionalWidth();
	//char finename[] = "data/two1.txt";
	//int totalNum = 500;
	//GenTwoDimensinoGridDataSet(finename, totalNum);

	char filename[] = "data/iris_4.txt";
	LoadData_line(filename);
	/*FILE *fp = fopen(filename, "r");
	char line[1024];
	while (fgets(line, 1024, fp))
	{
		puts(line);
	}
	fclose(fp);*/

	return 0;
}
