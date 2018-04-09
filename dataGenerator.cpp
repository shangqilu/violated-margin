#include "dataGenerator.h"
#include "headers.h"
#include "globalVar.h"
#include "tools.h"
using namespace std;

double GenUniformRandom(double rangeStart, double rangeEnd)
{
    if (rangeStart > rangeEnd) {
        PrintError("range start is larger than rangeEnd");
    }
	double r = rangeStart;
	r = r + (rangeEnd - rangeStart)*((double)rand()) / RAND_MAX;
	if (r < rangeStart) r = rangeStart;
	if (r > rangeEnd) r = rangeEnd;
    return r;
}


//using Box - Muller transform 
double GenGaussianRandom()
{
	double x1;
	do {
		x1 = GenUniformRandom(0.0, 1.0);
	} while (fabs(x1) < ZERO_ERROR);
	//because of log, x1 cannot be zero
	double x2 = GenUniformRandom(0.0, 1.0);
	double z = sqrt(-2.0 * log(x1)) * cos(2.0 * M_PI * x2);
	return z;
}

void GenRandomPoint(double rangeStart, double rangeEnd, Point &p)
{
	for (int i = 0; i < Dim; i++)
	{
		p.x[i] = GenUniformRandom(rangeStart, rangeEnd);
	}
}

void GenRandomPointInCircle(Point &center, double radius, Point &p)
{
	for (int i = 0; i < Dim; i++){
		double v = GenUniformRandom(-radius, radius);
		p.x[i] = center.x[i] + v;
	}
}

void GenRandomHyperPlane(double rangeStart, double rangeEnd, HyperPlane &plane)
{
	for (int i = 0; i < plane.d; i++){
		plane.w[i] = GenUniformRandom(rangeStart, rangeEnd);
	}
	plane.b = GenUniformRandom(rangeStart, rangeEnd);
}



void GenTwoDimensinoGridDataSet(char *filename, int totalNum)
{
	FILE *fp = fopen(filename, "w");
	if (fp == NULL) {
		PrintError("cannot open current file!");
	}

	int dimension = 2;
	double *x = new double[Dim];
	Point *cur = new Point(Dim, x);
	for (int i = 0; i < totalNum; i++)
	{
		
		GenRandomPoint(-1.0, 1.0, *cur);
		if (cur->x[0] * cur->x[1] > 0) {
			cur->y = 1;
		} else {
			cur->y = -1;
		}
		//write Point cur to file
		fprintf(fp, "%d", cur->y);
		for (int j = 0; j < dimension; j++)
		{
			fprintf(fp, " %d:%lf", j + 1, cur->x[j]);
		}
		fprintf(fp, "\n");
	}
	delete[]cur->x;
	delete cur;
	fclose(fp);
}

void GenMarginPoint(HyperPlane &plane, Point &center, double radius, double margin, Point &p)
{
	double *x = new double[Dim];
	Point cur(plane.d, x);
	while (1)
	{
		GenRandomPointInCircle(center, radius, cur);
		double value = Dot(plane.w, cur.x, plane.d) + plane.b;
		double denominator = 0;
		for (int i = 0; i < plane.d; i++) {
			denominator += plane.w[i] * plane.w[i];
		}
		double dis = value / sqrt(denominator);
		if (value > 0) {
			for (int i = 0; i < plane.d; i++) {
				p.x[i] = cur.x[i] + (margin - dis) * plane.w[i] / sqrt(denominator);
			}
			if (Distance(p, center, plane.d) < radius) {
				p.y = 1;
				break;
			}
			
		}
		else if (value < 0){
			for (int i = 0; i < plane.d; i++) {
				p.x[i] = cur.x[i] - (margin - dis) * plane.w[i] / sqrt(denominator);
			}
			if (Distance(p, center, Dim) < radius) {
				p.y = -1;
				break;
			}
		}
		
	}
	delete[]cur.x;
}

void GenMarginDataSet(char *filename, double margin, double radius,
	int totalNum, int noiseNum)
{
	cout << "totalNum: " << totalNum << " noiseNum:" << noiseNum << endl;
	cout << "radius: " << radius << " margin:" << margin << endl;
	if (margin > radius) {
		PrintError("Error parameters: margin is greater than radius");
	}
	FILE *fp = fopen(filename, "w");
	if (fp == NULL) {
		PrintError("cannot open current file!");
	}
	HyperPlane plane(Dim);
	//generate a hyperplane with wights and bias ranging from 0 to 1
	GenRandomHyperPlane(0.0, 1.0, plane);
	int len = strlen(filename);
	char *planefile = new char[len+10];
	strcpy(planefile, filename);
	strcpy(planefile + len, ".discribe");
	planefile[len + 9] = '\0';
	FILE *f_plane = fopen(planefile, "w");
	if (f_plane == NULL) {
		PrintError("can not open the plane file");
	}
	fprintf(f_plane, "optimal plane:");
	for (int j = 0; j < Dim; j++)
	{
		fprintf(f_plane, " %lf", plane.w[j]);
	}
	fprintf(f_plane, " b: %lf\n", plane.b);
	fprintf(f_plane, "radius: %lf\n", radius);
	fprintf(f_plane, "margin: %lf\n", margin);
	fprintf(f_plane, "totalNum: %d\n", totalNum);
	fprintf(f_plane, "noiseNum: %d\n", noiseNum);
	delete[] planefile;

	//choos a point in the plane as a center point
	//choose one axis, and set the coordinates in another dimension with value 1
	double *x = new double[Dim];
	Point center(Dim, x);
	for (int i = Dim - 1; i >= 0; i--) {
		double *x = new double[Dim];
		Point axis(Dim, x);
		axis.x[i] = 1;
		if (Dot(plane.w, axis.x, Dim) > ZERO_ERROR) {
			double sum = 0;
			for (int j = 0; j < Dim; j++) {
				if (j != i) {
					center.x[j] = 1;
					sum += plane.w[j];
				}
			}
			sum += plane.b;
			center.x[i] = -sum / plane.w[i];
			delete[]axis.x;
			break;
		}
		delete[]axis.x;
	}
	puts("centering at:");
	for (int i = 0; i < Dim; i++)
	{
		cout << center.x[i] << " ";
	}
	cout << endl;
	fprintf(f_plane, "center:");
	for (int j = 0; j < Dim; j++)
	{
		fprintf(f_plane, " %lf", center.x[j]);
	}
	fprintf(f_plane, "\n");
	fclose(f_plane);

	//generating points within a ball
	int rightNum = totalNum - noiseNum;
	//generating O(d) points at the margin 
	int marginNum = Dim * 5;
	int times = rightNum / marginNum;
	int cnt = 0;
	x = new double[Dim];
	Point cur(Dim, x);
	for (int i = 0; i < rightNum; i++) {
		if (i % times == 0)
		{
			GenMarginPoint(plane, center, radius, margin, cur);
			cnt++;
		}
		else{
			while (1)
			{
				GenRandomPointInCircle(center, radius, cur);
				if (Distance(plane, cur, Dim) > margin) {
					double value = Dot(plane.w, cur.x, Dim) + plane.b;
					if (value > 0) {
						cur.y = 1;
					}
					else {
						cur.y = -1;
					}
					break;
				}
			}
		}
		//write Point cur to file
		fprintf(fp, "%d", cur.y);
		for (int j = 0; j < Dim; j++)
		{
			fprintf(fp, " %d:%lf", j + 1, cur.x[j]);
		}
		fprintf(fp, "\n");
	}
	printf("generating %d margin points\n", cnt);

	//generating noise points
	for (int i = 0; i < noiseNum; i++) {
		while (1)
		{
			GenRandomPointInCircle(center, radius, cur);
			if (Distance(plane, cur, plane.d) > margin) {
				double value = Dot(plane.w, cur.x, Dim) + plane.b;
				if (value > 0) {
					cur.y = -1;
				}
				else {
					cur.y = 1;
				}
				break;
			}
		}
		//write Point cur to file
		fprintf(fp, "%d", cur.y);
		for (int j = 0; j < Dim; j++)
		{
			fprintf(fp, " %d:%lf", j + 1, cur.x[j]);
		}
		fprintf(fp, "\n");
	}

	delete[]center.x;
	delete[]cur.x;

	fclose(fp);
	return;
}