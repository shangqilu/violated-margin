#ifndef __SIMPLEX_H__

#define __SIMPLEX_H__
#include "headers.h"

using namespace std;

/*
*   a data structure to store the linear program in standard form
*   the objective function is to maximize C*x
*   subject to Ax <= b
*/
struct Simplex_Node{
    double **A; //store the linear program in standard form
    double *b;  //Ax <= b
    double *C;  //the objective function vector

    int *N;     //non-basic variables
    int *B;     //basic variables

    int m;      // number of constraints
    int n;      // number of features
    double v;   //objective function value

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
    
    ~Simplex_Node()
    {
        if (A != NULL) {
            for (int i = 0; i < m; i++)
            {
                delete []A[i];
            }
            delete []A;
            A = NULL;
        }
        if (b != NULL) {
            delete []b;
            b = NULL;
        }
        if (C != NULL) {
            delete []C;
            C = NULL;
        }
        if (N != NULL) {
            delete []N;
            N = NULL;
        }
        if (B != NULL) {
            delete []B;
            B = NULL;
        }
    }
    
};

/*
*   a structure to store the result of simplex algorithm
*/
struct LPresult{
    double *x;      //the solution to a LP
    int flag;       //flag = 1 find an optimal solution
                    //flag = 0 infeasible or unbounded
    double value;   //optimal objective function value
    LPresult(int d)
    {
        this->x = new double[d];
        for (int i = 0; i < d; i++) this->x[i] = 0;
        this->flag = 0;
        this->value = 0;
    }

    ~LPresult()
    {
        if (x != NULL) {
            delete []x;
            x = NULL;
        }
    }
};

/*
*   compute the solution to a standard linear program
*/
LPresult Simplex(double **input_A, double *input_b, double *input_C, int m, int n);
/*
*   pivot operation: given an entering variable and a leaving variable
*   exchange their roles
*/
void Pivot(Simplex_Node &node, int leaving, int entering);

/*
*   determine whether current instance is feasible
*   if feasible return true
*       return a Simplex node with the initial basic solution feasible
*   else return false
*/
bool Initial_Simplex(Simplex_Node &node);

void PrintLPresult(LPresult result, int dimension);
void PrintSimplexNode(Simplex_Node node);

/*
*   find the hyperplane with largest gap
*   along the direction of the i th dimension coordinate axis
*/
bool OneDirectionLPClassification(PointSet trainPoints, HyperPlane &plane, int dimension, int direction);

/*
*   compute the 1/sqrt(dimension)-approximation hyperplane through liner programming
*/
bool LPclassification(PointSet trainPoints, HyperPlane &plane, int dimension);
#endif // __SIMPLEX_H__
