#include "simplex.h"
#include "headers.h"
#include "globalVar.h"
#include "tools.h"

LPresult Simplex(double **input_A, double *input_b, double *input_C, int m, int n)
{
    //puts("RunSimplex...");
    Simplex_Node node(input_A, input_b, input_C, m, n);
    LPresult result = LPresult(n);
	int num;
    bool re = Initial_Simplex(node, num);
	/*if (num == 100000) {
		char filename[] = "data/cycleLP.txt";
		FILE *fp = fopen(filename, "w");
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++)
			{
				fprintf(fp, "%lf ", input_A[i][j]);
			}
			fprintf(fp, "\n");
		}
		for (int i = 0; i < m; i++) fprintf(fp, "%lf ", input_b[i]);
		fprintf(fp, "\n");
		for (int i = 0; i < n; i++) fprintf(fp, "%lf ", input_C[i]);
		fprintf(fp, "\n");
		exit(0);
	}*/
	//puts("after initial");
    if (!re) {
        //puts("It is infeasible in initial simplex!");
        return result;
    }

    while(1){
        //choose the variable with the smallest index
        //puts("pivoting...");
        //PrintSimplexNode(node);
        int min_index = n + m; //non basic variable starting form m
        int entering = n + m;
        for (int j = 0; j < n; j++) {
			if ((node.C[j] - 0)> ZERO) {
                //puts("can not stop");
                //if (min_index > node.N[j]) {
                    //cout << "asdf" << endl;
                    min_index = node.N[j];
                    entering = j;
				//}
					break;
            }
        }
        //printf("choosing entering: %d %d\n", entering, min_index);
        if(entering == n + m) {
            //end computing with an optimal solution
            //puts("find the optimal value");
            result.flag = 1;
            break;
        }
        double max_delta = MAX_DOUBLE;
        min_index = m + n;
        int leaving = m + n;
        for (int i = 0; i < m; i++) {
			if (node.A[i][entering] - 0 > ZERO_ERROR) {
                double tmp = node.b[i]/node.A[i][entering];
                //printf("*** %d %lf\n", i, tmp);
                if (max_delta > tmp) {
                    //printf("update leaving %d %lf %lf", i, tmp);
                    leaving = i;
                    max_delta = tmp;
                    min_index = node.B[i];
                }
				/*if (fabs(max_delta - tmp)< ZERO_ERROR) {
                    if(min_index > node.B[i]) {
                        leaving = i;
                        max_delta = tmp;
                        min_index = node.B[i];
                    }
                }*/
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
    tmp = node.C[entering]*node.b[leaving];
    node.v += tmp;
    for (int j = 0; j < node.n; j++) {
        if (j != entering) {
            double tmp = node.C[entering]*node.A[leaving][j];
            node.C[j] -= tmp;
            //cout << tmp << " " << node.C[j] << endl;
        }
    }
    tmp = node.C[entering]*node.A[leaving][entering];
    node.C[entering] = - tmp;
    //compute the new basic variables
    int tmp_entering_index = node.N[entering];
    node.N[entering] = node.B[leaving];
    node.B[leaving] = tmp_entering_index;
}


bool Initial_Simplex(Simplex_Node &node, int &num)
{
    //puts("Run Initial_Simplex...");
    int k = -1;
    double min_value = MAX_DOUBLE;

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

    double *back_c = new double[n-1];
    //change the objective function into -x_n
    for (int j = 0; j < n-1; j++) {
        back_c[j] = node.C[j];
        node.C[j] = 0;
    }
    node.C[n-1] = -1;
    //PrintSimplexNode(node);
    int leaving = k;
    Pivot(node, leaving, n-1);
    //PrintSimplexNode(node);
    //L_aux find the optimal solution of L_aux
    int flag = 0;
	num = 0;
    while(1){
        //choose the variable with the smallest index
		num++;
		/*if (num == 100000) {
			delete []back_c;
			return false;
		}*/
		//printf("Initial pivoting... %lf\n", node.v);
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
					//break;
            }
        }
        //printf("choosing entering: %d %d\n", entering, min_index);
        if(entering == n + m) {
            //end computing with an optimal solution
            //puts("find the optimal value in Initial-simplex");
            //PrintSimplexNode(node);
            flag = 1;
            break;
        }
        double max_delta = MAX_DOUBLE;
        min_index = m + n;
        int leaving = m + n;
        for (int i = 0; i < m; i++) {
            if (node.A[i][entering]-0 > ZERO) {
                double tmp = node.b[i]/node.A[i][entering];
                //printf("*** %d %lf\n", i, tmp);
                if (max_delta > tmp) {
                    //printf("update leaving %d %lf %lf", i, tmp);
                    leaving = i;
                    max_delta = tmp;
                    min_index = node.B[i];
                }
				if (fabs(max_delta - tmp)< ZERO_ERROR) {
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
            puts("It is unbounded in LP_aux!");
            break;
        } else{
            Pivot(node, leaving, entering);
        }

    }
	if (!flag || fabs(node.v - 0) > ZERO_ERROR) {
        //cout << flag << " " << node.v << endl;
        //puts("infeasible!");
        delete []back_c;

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
    //for (int i = 0; i < n; i++)
    //    cout << back_c[i] << " ";
    //cout << endl;
    //puts("*****************");
    delete back_c;
    return true;
}

void PrintLPresult(LPresult &result, int dimension)
{
    printf("objective value: %lf\n", result.value);
    printf("X: ");
    for (int i = 0; i < dimension; i++)
    {
        printf("%lf ", result.x[i]);
    }
    cout << endl;
}

void PrintSimplexNode(Simplex_Node &node)
{
    int m = node.m;
    int n = node.n;
    puts("Matrix A:");
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << node.A[i][j] << " ";
        }
        cout << endl;
    }
    puts("b: ");
    for (int i = 0; i < m; i++) {
        cout << node.b[i] << " ";
    }
    cout << endl;
    puts("objective function C:");
    for (int j = 0; j < n; j++) {
        cout << node.C[j] << " ";
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
    printf("Cur objecrive value :");
    cout << node.v << endl;
}

//the direction range from 1 to d means d directions along the axes
//if direction = 0, there are no direction constraints
bool OneDirectionLPClassification(PointSet &trainPoints, HyperPlane &plane, int direction)
{
    int number_constraints, number_variables;
	number_variables = 2 * (Dim + 1) + 1;
    if (direction == 0)   //simple LP, no specific direction
    {
        number_constraints = trainPoints.size();
    }
    else
    {
		number_constraints = trainPoints.size() + 2 * Dim;
    }
    double **input_A, *input_b, *input_C;
    input_A = new double*[number_constraints];
    for (int i = 0; i < number_constraints; i++)
    {
        input_A[i] = new double[number_variables];
    }
    input_b = new double[number_constraints];
    input_C = new double[number_variables];
    int n = trainPoints.size();

    for (int i = 0;i < number_constraints; i++)
    {
        for (int j = 0; j < number_variables; j++) {
            input_A[i][j] = 0;
        }
    }
    for (int i = 0; i < n; i++)
    {
        int cur_label = trainPoints[i]->y;//label: 1 or -1
		for (int j = 0; j < Dim; j++)
        {
            input_A[i][2*j] = - trainPoints[i]->x[j] * cur_label;
            input_A[i][2*j+1] = trainPoints[i]->x[j] * cur_label;
        }
		input_A[i][2 * Dim] = -1 * cur_label;
		input_A[i][2 * Dim + 1] = 1 * cur_label;
		input_A[i][2 * Dim + 2] = 1;//always 1 in last dimension
        input_b[i] = 0;
    }
    int cur_row = n;
    if(direction > 0)   //add constraints information about w
    {
		for (int j = 0; j < Dim; j++) {
            if (j == direction-1) {
                input_A[cur_row][2*j] = 1;// w_j<=1
                input_A[cur_row][2*j+1] = -1;
                input_b[cur_row] = 1;

                input_A[cur_row+1][2*j] = -1; //w_j>=1
                input_A[cur_row+1][2*j+1] = 1;
                input_b[cur_row+1] = -1;
            }
            else {
                input_A[cur_row][2*j] = 1;// w_j<=1
                input_A[cur_row][2*j+1] = -1;
                input_b[cur_row] = 1;

                input_A[cur_row+1][2*j] = -1; //w_j>=-1
                input_A[cur_row+1][2*j+1] = 1;
                input_b[cur_row+1] = 1;
            }
            cur_row += 2;
        }
    }
    for (int j = 0; j < number_variables - 1; j++)
    {
        input_C[j] = 0;
    }
    input_C[number_variables-1] = 1;
    
    LPresult result = Simplex(input_A, input_b, input_C, number_constraints, number_variables);

    for (int i = 0; i < number_constraints; i++) {
        delete []input_A[i];
    }
    delete []input_A;
    delete []input_b;
    delete []input_C;


    if (result.flag == 0) {
        return false;
    }
    //printf("number of constraints: %d\n", number_constraints);
    //PrintLPresult(result, number_variables);
	for (int j = 0; j < Dim; j++) {
        plane.w[j] = result.x[2*j]-result.x[2*j+1];
    }
	plane.b = result.x[2 * Dim] - result.x[2 * Dim + 1];
    return true;
}


bool LPclassification(PointSet &trainPoints, HyperPlane &plane)
{
    HyperPlane cur_plane(Dim);
    double max_margin = 0;
    int flag = 0;
	for (int j = 0; j < Dim; j++) {
		bool found = OneDirectionLPClassification(trainPoints, cur_plane, j + 1);
        if (found) {
            flag = 1;
            double cur_margin = 0;
            bool separable = MinimumSeparableDistance(trainPoints, cur_plane, cur_margin);
            if (!separable) {
                //printf("wrong in simplex, current dimension %d\n", j);
                //minimumseparabledistance has a more serious condition to check whether its separable.
                //PrintHyperPlane(cur_plane, dimension);
                flag = 0;

            }
            if (separable && cur_margin > max_margin) {
                //printf("find a solution in dimension %d with margin %lf\n", j, cur_margin);
                //PrintHyperPlane(cur_plane);
                max_margin = cur_margin;
                
				for (int j = 0; j < Dim; j++) {
                    plane.w[j] = cur_plane.w[j];
                }
                plane.b = cur_plane.b;
                //CopyHyperPlane(plane, cur_plane);
            }
        }

    }
    if (!flag) return false;
    else return true;
}
