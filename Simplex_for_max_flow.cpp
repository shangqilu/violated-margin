#include <cstring>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <map>
using namespace std;

#define MAX_DOUBLE 1e9
#define ZERO  1e-6

struct Simplex_Node{
    int **A;
    int *b;
    int *C;

    int *N; //non-basic variables
    int *B; //basic variables

    int m;// number of constraints
    int n;// number of features
    int v; //objective function value

    Simplex_Node(int **input_A, int *input_b, int *input_C, int m, int n)
    {
        this->A = new int*[m];
        for (int i = 0; i < m; i++) {
            this->A[i] = new int[n+1];
        }
        this->b = new int[m];
        this->C = new int[n+1];
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
    int *x;
    int flag;
    int value;
    LPresult(int d)
    {
        this->x = new int[d];
        this->flag = 0;
        this->value = 0;
    }
};


LPresult Simplex(int **input_A, int *input_b, int *input_C, int m, int n);
void Pivot(Simplex_Node &node, int leaving, int entering);
bool Initial_Simplex(Simplex_Node &node);
void PrintLPresult(LPresult result, int dimension);
void PrintSimplexNode(Simplex_Node node);


LPresult Simplex(int **input_A, int *input_b, int *input_C, int m, int n)
{
    //puts("RunSimplex...");
    Simplex_Node node = Simplex_Node(input_A, input_b, input_C, m, n);
    LPresult result = LPresult(n);
    bool re = Initial_Simplex(node);
    if (!re) {
        //puts("It is infeasible!");
        return result;
    }

    while(1){
        //choose the variable with the smallest index
        //puts("pivoting...");
        //PrintSimplexNode(node);
        int min_index = n + m; //non basic variable starting form m
        int entering = n + m;
        for (int j = 0; j < n; j++) {
            if((node.C[j]-0)> ZERO) {
                //puts("can not stop");
                if (min_index > node.N[j]) {
                    //cout << "asdf" << endl;
                    min_index = node.N[j];
                    entering = j;
                }
            }
        }
        //printf("choosing entering: %d %d\n", entering, min_index);
        if(entering == n + m) {
            //end computing with an optimal solution
            //puts("find the optimal value");
            result.flag = 1;
            break;
        }
        int max_delta = MAX_DOUBLE;
        min_index = m + n;
        int leaving = m + n;
        for (int i = 0; i < m; i++) {
            if (node.A[i][entering] > 0) {
                int tmp = node.b[i]/node.A[i][entering];
                //printf("*** %d %d\n", i, tmp);
                if (max_delta > tmp) {
                    //printf("update leaving %d %d %d", i, tmp);
                    leaving = i;
                    max_delta = tmp;
                    min_index = node.B[i];
                }
                if (fabs(max_delta - tmp)< 1e-6) {
                    if(min_index > node.B[i]) {
                        leaving = i;
                        max_delta = tmp;
                        min_index = node.B[i];
                    }
                }
            }
        }
        //printf("choosing leaving: %d %d\n", leaving, min_index);
        if(leaving == m + n) {
            //unbounded
            //puts("It is unbounded!");
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
    //compute new coefficients for the new equation for new basic variable x_e
    int tmp = node.A[leaving][entering];
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
    int tmp_entering_index = node.N[entering];
    node.N[entering] = node.B[leaving];
    node.B[leaving] = tmp_entering_index;
}


bool Initial_Simplex(Simplex_Node &node)
{
    //puts("Run Initial_Simplex...");
    int k = -1;
    int min_value = MAX_DOUBLE;

    //printf("m: %d n: %d\n", node.m, node.n);
    for (int i = 0;i < node.m; i++) {
        if(min_value > node.b[i])
        {
            min_value = node.b[i];
            k = i;
        }
    }
    //The initial basic solution is feasible
    if (min_value >= 0) {
        //puts("initial basic solution ok ");
        return true;
    }
    //add a new non basic variable x_n
    node.n++;
    int m = node.m;
    int n = node.n;
    //puts("constructing LP_aux ...");
    for (int i = 0; i < m; i++) {
        node.A[i][n-1] = -1;
    }
    node.N[n-1] = n-1+m;

    int *back_c = new int[n-1];
    //change the objective function into -x_n
    for (int j = 0; j < n-1; j++) {
        back_c[j] = node.C[j];
        node.C[j] = 0;
    }
    node.C[n-1] = -1;
    PrintSimplexNode(node);
    int leaving = k;
    Pivot(node, leaving, n-1);

    //L_aux find the optimal solution of L_aux
    int flag = 0;
    while(1){
        //choose the variable with the smallest index
        //puts("pivoting...");
        //PrintSimplexNode(node);
        int min_index = n + m; //non basic variable starting form m
        int entering = n + m;
        for (int j = 0; j < n; j++) {
            if((node.C[j]-0)> ZERO) {
                //puts("can not stop");
                if (min_index > node.N[j]) {
                    //cout << "asdf" << endl;
                    min_index = node.N[j];
                    entering = j;
                }
            }
        }
        //printf("choosing entering: %d %d\n", entering, min_index);
        if(entering == n + m) {
            //end computing with an optimal solution
            //puts("find the optimal value");
            flag = 1;
            break;
        }
        int max_delta = MAX_DOUBLE;
        min_index = m + n;
        int leaving = m + n;
        for (int i = 0; i < m; i++) {
            if (node.A[i][entering] > 0) {
                int tmp = node.b[i]/node.A[i][entering];
                //printf("*** %d %d\n", i, tmp);
                if (max_delta > tmp) {
                    //printf("update leaving %d %d %d", i, tmp);
                    leaving = i;
                    max_delta = tmp;
                    min_index = node.B[i];
                }
                if (fabs(max_delta - tmp)< 1e-6) {
                    if(min_index > node.B[i]) {
                        leaving = i;
                        max_delta = tmp;
                        min_index = node.B[i];
                    }
                }
            }
        }
        //printf("choosing leaving: %d %d\n", leaving, min_index);
        if(leaving == m + n) {
            //unbounded
            //puts("It is unbounded!");
            break;
        } else{
            Pivot(node, leaving, entering);
        }

    }
    if (!flag || fabs(node.v - 0) > ZERO) {
        //puts("infeasible!");
        return false;
    }
    //if x_n is a basic variable, we should do a pivot
    //find any a_n_e != 0;
    for (int i = 0; i < m; i++) {
        if (node.B[i] == n-1+m) {
            for (int j = 0; j < n; j++) {
                if(node.A[i][j] != 0) {
                    Pivot(node, i, j);
                    break;
                }
            }
            break;
        }
    }
    //swap x_n with the last column and remove x_n
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

    //puts("*****************");
    //recover the objective function
    for (int j = 0; j < n; j++)
    {
        if (node.N[j] >= m) {
            node.C[j] = back_c[node.N[j]-m];
        }else {
            node.C[j] = 0;
        }

    }
    //PrintSimplexNode(node);
    //replace the basic variables in the objective function
    for (int i = 0; i < m; i++) {
        if (node.B[i] >= m) {
            for (int j = 0; j < n; j++) {
                node.C[j] += -1*node.A[i][j]*back_c[node.B[i]-m];
            }
            node.v += node.b[i]*back_c[node.B[i]-m];
        }
    }
    for (int i = 0; i < n; i++)
        cout << back_c[i] << " ";
    cout << endl;
    //puts("*****************");
    delete back_c;
    return true;
}

void PrintLPresult(LPresult result, int dimension)
{
    printf("objective value: %d\n", result.value);
    printf("X: ");
    for (int i = 0; i < dimension; i++)
    {
        printf("%d ", result.x[i]);
    }
    cout << endl;
}

void PrintSimplexNode(Simplex_Node node)
{
    int m = node.m;
    int n = node.n;
    puts("Matrix A:");
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d ", node.A[i][j]);
        }
        cout << endl;
    }
    puts("b: ");
    for (int i = 0; i < m; i++) {
        printf("%d ", node.b[i]);
    }
    cout << endl;
    puts("objective function C:");
    for (int j = 0; j < n; j++) {
        printf("%d ", node.C[j]);
    }
    cout << endl;
    puts("Basic variables:");
    for (int i = 0; i < m; i++) {
        printf("%d ", node.B[i]);
    }
    cout << endl;
    puts("Non basic variables:");
    for (int j = 0; j < n; j++) {
        printf("%d ", node.N[j]);
    }
    cout << endl;
    printf("Cur objecrive value %d\n", node.v);
}

const int mmax = 202;
int Map[mmax][mmax];
int A[mmax*3][mmax*2];
int input_B[mmax*3];
int input_C[2*mmax];
int *input_A[mmax*3];

map<int, int> index_map;
int main()
{
    freopen("data/LPtest/input.txt", "r", stdin);
    int m, n;
    while(cin >> n >> m) {
        if (n == 0) {
            puts("0");
            continue;
        }
        index_map.clear();
        memset(Map, 0, sizeof(Map));

        for (int i = 0; i < mmax*3; i++) {
            input_A[i] = A[i];
        }


        for (int i = 0; i < mmax*3; i++) {
            for (int j = 0;j < 2*mmax; j++)
                input_A[i][j] = 0;
        }


        for (int i = 0; i < mmax*3; i++) {
            input_B[i] = 0;
        }
        for (int j = 0; j < 2*mmax; j++) {
            input_C[j] = 0;
        }
        int u, v, cost;
        for (int i = 0; i < n; i++) {
            cin >> u >> v >> cost;
            Map[u-1][v-1] += cost;
        }

        int column = 0;
        int cnt = 0;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                if (Map[i][j] > 0) {
                    int index = i*m+j;
                    if (!index_map.count(index)) {
                        index_map[index] = column++;
                    }
                    input_A[cnt][index_map[index]] = 1;
                    input_B[cnt++] = Map[i][j];
                }
            }
        }

        for (int i = 1; i < m-1; i++) {
            for (int j = 0; j < m; j++) {
                if (j!=i) {
                    if(Map[i][j] > 0) {
                        if (!index_map.count((i*m+j)))index_map[i*m+j] = column++;
                        input_A[cnt][index_map[i*m+j]] = 1;
                        input_A[cnt+1][index_map[i*m+j]] = -1;
                    }
                    if (Map[j][i] > 0) {
                        if (!index_map.count((j*m+i)))index_map[j*m+i] = column++;
                        input_A[cnt][index_map[j*m+i]] = -1;
                        input_A[cnt+1][index_map[j*m+i]] = 1;
                    }
                }
            }
            input_B[cnt] = input_B[cnt+1] = 0;
            cnt += 2;
        }
        for (int i = 0;i < m; i++)
        {
           if (Map[0][i] > 0) input_C[index_map[i]] = 1;
        }
        LPresult result = Simplex(input_A, input_B, input_C, cnt, m*m);
        printf("%d\n", result.value);
        //PrintLPresult(result, m*m);
    }
    return 0;
}
