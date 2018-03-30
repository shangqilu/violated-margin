#ifndef __SIMPLEX_H__

#define __SIMPLEX_H__
#include "headers.h"

using namespace std;

struct Simplex_Node{
    double **A;
    double *b;
    double *C;

    int *N; //non-basic variables
    int *B; //basic variables

    int m;// number of constraints
    int n;// number of features
    double v; //objective function value

    Simplex_Node(double **input_A, double *input_b, double *input_C, int m, int n)
    {
        this->A = new double*[m];
        for (int i = 0; i < m; i++) {
            this->A[i] = new double[n+1];
        }
        this->b = new double[m];
        this->C = new double[n+1];
        this->N = new int[n+1];
        this->B = new int[m];

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                this->A[i][j] = input_A[i][j];
            }
            this->b[i] = input_b[i];
            this->B[i] = i;
        }
        for (int j = 0; j < n; j++) {
            this->C[j] = input_C[j];
            this->N[j] = j+m;
        }
        this->v = 0;
        this->m = m;
        this->n = n;
    }
};

struct LPresult{
    double *x;
    int flag;
    double value;
    LPresult(int d)
    {
        this->x = new double[d];
        this->flag = 0;
        this->value = 0;
    }
};


LPresult Simplex(double **input_A, double *input_b, double *input_C, int m, int n);
void Pivot(Simplex_Node &node, int leaving, int entering);
bool Initial_Simplex(Simplex_Node &node);
void PrintLPresult(LPresult result, int dimension);
void PrintSimplexNode(Simplex_Node node);
bool OneDirectionLPClassification(PointSet trainPoints, HyperPlane &plane, int dimension, int direction);
bool LPclassification(PointSet trainPoints, HyperPlane &plane, int dimension);
#endif // __SIMPLEX_H__
