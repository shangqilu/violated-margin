#include "simplex.h"
#include "headers.h"

LPresult Simplex(double **input_A, double *input_b, double *input_C, int m, int n)
{
    Simplex_Node node = Simplex_Node(input_A, input_b, input_C, m, n);
    LPresult result = LPresult(n);
    bool re = Initial_Simplex(node);
    if (!re) {
        puts("It is infeasible!");
        return result;
    }

    while(1){
        //choose the variable with the smallest index
        int min_index = n;
        int entering = n;
        for (int j = 0; j < n; j++) {
            if(node.C[j] > 0) {
                if (min_index > node.N[j]) {
                    min_index = node.N[j];
                    entering = j;
                }
            }
        }
        if(entering == n) {
            //end computing with an optimal solution
            result.flag = 1;
            break;
        }
        double max_delta = MAX_DOUBLE;
        min_index = m;
        int leaving = m;
        for (int i = 0; i < m; i++) {
            if (node.A[i][entering] > 0) {
                double tmp = node.b[i]/node.A[i][entering];
                if (max_delta > tmp) {
                    leaving = i;
                    min_index = node.B[i];
                }else if (max_delta == tmp) {
                    if(min_index > node.B[i]) {
                        leaving = i;
                        min_index = node.B[i];
                    }
                }
            }
        }
        if(leaving == m) {
            //unbounded
            puts("It is unbounded!");
            return result;
        } else{
            Pivot(node, leaving, entering);
        }

    }
    //find the solution
    for (int i = 0; i < m; i++) {
        if (node.B[i] >= m) {
            result.x[node.B[i]-m] = node.b[i];
        }
    }
    result.value = node.v;
    return result;
}




void Pivot(Simplex_Node &node, int leaving, int entering)
{
    double **new_A;
    //compute new coefficients for the new equation for new basic variable x_e
    double tmp = node.A[leaving][entering];
    node.b[leaving] /= tmp;
    for (int j = 0; j < node.n; j++)
    {
        if(j != entering) {
            node.A[leaving][j] /= tmp;
        }
    }
    node.A[leaving][entering] = 1/tmp;
    //compute remaining constraints
    for (int i = 0; i < node.m; i++) {
        if(i != leaving) {
            node.b[i] -= node.A[i][entering]* node.b[leaving];
            for (int j = 0; j < node.n; j++) {
                if (j != entering) {
                    node.A[i][j] -= node.A[i][entering] * node.A[leaving][j];
                }
            }
            node.A[i][entering] = - node.A[i][entering] * node.A[leaving][entering];
        }
    }
    //compute the objective function
    node.v += node.C[entering]*node.b[leaving];
    for (int j = 0; j < node.n; j++) {
        if (j != entering) {
            node.C[j] -= node.C[entering]*node.A[leaving][j];
        }
    }
    node.C[entering] = -node.C[entering]*node.A[leaving][entering];
    //compute the new basic variables
    node.N[entering] = leaving;
    node.B[leaving] = entering;
}


bool Initial_Simplex(Simplex_Node &node)
{
    int k = -1;
    double min_value = MAX_DOUBLE;
    for (int i = 0;i < node.m; i++) {
        if(min_value > node.b[i])
        {
            min_value = node.b[i];
            k = i;
        }
    }
    //The initial basic solution is feasible
    if (min_value >= 0) {
        return true;
    }
    //add a new non basic variable x_n
    node.n++;
    int m = node.m;
    int n = node.n;

    for (int i = 0; i < m; i++) {
        node.A[i][n-1] = -1;
    }
    node.N[n-1] = n-1+m;

    double *back_c = new double[n-1];
    //change the objective function into -x_n
    for (int j = 0; j < n-1; j++) {
        back_c[j] = node.C[j];
        node.C[j] = 0;
    }
    node.C[n-1] = -1;

    int leaving = k;

    Pivot(node, leaving, n-1);
    //L_aux find the optimal solution of L_aux
    int flag = 0;
    while(1){
        //choose the variable with the smallest index
        int min_index = n;
        int entering = n;
        for (int j = 0; j < n; j++) {
            if(node.C[j] > 0) {
                if (min_index > node.N[j]) {
                    min_index = node.N[j];
                    entering = j;
                }
            }
        }
        if(entering == n) {
            //end computing with an optimal solution
            flag = 1;
            break;
        }
        double max_delta = MAX_DOUBLE;
        min_index = m;
        leaving = m;
        for (int i = 0; i < m; i++) {
            if (node.A[i][entering] > 0) {
                double tmp = node.b[i]/node.A[i][entering];
                if (max_delta > tmp) {
                    leaving = i;
                    min_index = node.B[i];
                }else if (max_delta == tmp) {
                    if(min_index > node.B[i]) {
                        leaving = i;
                        min_index = node.B[i];
                    }
                }
            }
        }
        if(leaving == m) {
            //unbounded
            puts("It is unbounded!");
            break;
        } else{
            Pivot(node, leaving, entering);
        }

    }
    if (!flag || node.v != 0) {
        puts("infeasible!");
        return false;
    }
    //if x_n is a basic variable
    for (int i = 0; i < m; i++) {
        if (node.B[i] == n-1+m) {
            for (int j = 0; j < n; j++) {
                if(node.A[i][j] != 0) {
                    Pivot(node, i, j);
                }
            }
            break;
        }
    }
    for (int j = 0; j < n; j++) {
        if (node.N[j] == n-1+m && j != n-1) {
            for (int i = 0; i < m; i++) {
                node.A[i][j] = node.A[i][n-1];
            }
            node.N[j] = node.N[n-1];

        }
    }
    n--;
    node.n--;
    for (int j = 0; j < n; j++)
    {
        if (node.N[j] >= m) {
            node.C[j] = back_c[node.N[j]-m];
        }else {
            node.C[j] = 0;
        }

    }
    //replace the basic variables in the objective function
    for (int i = 0; i < m; i++) {
        if (node.N[i] >= m) {
            for (int j = 0; j < n; j++) {
                node.C[j] += node.A[i][j]*back_c[node.N[i]-m];
            }
            node.v += node.b[i]*back_c[node.N[i]-m];
        }
    }
    return true;
}
