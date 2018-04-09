#include "directionalWidth.h"
#include "globalVar.h"
#include "tools.h"

using namespace std;

bool TwoApproxiDiameter(PointSet &points, int dimension, int &s, int &t)
{
	int n = points.size();
	int cnt = 0;
	double *direction = new double[dimension];
	double *axis = new double[dimension];
	for (int i = 0; i < dimension - 1; i++) axis[i] = 0;
	axis[dimension - 1] = 1;
	while (1)
	{
		s = rand() % n;
		double max_dis = -MAX_DOUBLE;
		t = -1;
		for (int i = 0; i < n; i++)
		{
			if (i != s)
			{
				double cur_dis = Distance(*points[s], *points[i], dimension);
				if (cur_dis > max_dis)
				{
					t = i;
					max_dis = cur_dis;
				}
			}
		}
		int flag = 0;
		//Point direction = PointMinus(points[s], points[t], dimension);
		
		for (int i = 0; i < dimension; i++) direction[i] = points[s]->x[i] - points[t]->x[i];
		//the hyper plane s and t defined is not parallel to last dimension coordinate axis
		if (fabs(Dot(direction, axis, dimension)) > ZERO_ERROR)
		{
			break;
		}
		cnt++;
		if (cnt == MaxIterTwoApprox)
		{
			break;
		}
	}
	delete[]direction;
	delete[]axis;
	if (cnt == MaxIterTwoApprox)
	{
		return false;
	}
	else
		return true;
}




bool RecursionMinimumBoudingBox(PointSet &points, BoundingBox &curbox, double **MainTainT, int dimension)
{
	if (dimension == 0)
	{
		return true;
	}
	//puts("Computing the approximate minimum bouding box recursively");
	//printf("cur dimension %d...\n", dimension);

	int n = points.size();
	if (n < 2)
	{
		//puts("there are no more than 2 points!");
		return false;
	}
	int s, t;
	double *pt_x = new double[dimension];
	Point d_axis(dimension, pt_x);
	bool found = TwoApproxiDiameter(points, dimension, s, t);
	if (!found)
	{
		//puts("points are in the same coordinates in current dimension!!! ");
		//puts("data can be presented in lower dimension");
		//puts("choose d dimension axis as the normal vector");
		//Point axis = Point(dimension);   //problem !!!
		//axis.x[dimension-1] = 1;
		d_axis.x[dimension - 1] = 1;

	}
	else
	{
		//d_axis = PointMinus(points[s], points[t], dimension);
		for (int i = 0; i < dimension; i++) d_axis.x[i] = points[s]->x[i] - points[t]->x[i];
	}
	//printf("current s and t: %d %d\n", s, t);
	//all points is above point t

	//find those d-1 basis

	double **CurT = new double*[Dim + 1];
	//store those d basis vectors
	for (int i = 0; i < Dim + 1; i++)
	{
		CurT[i] = new double[Dim + 1];
		for (int j = 0; j < Dim + 1; j++)
		{
			CurT[i][j] = 0;
		}
	}
	//filling the transforming matrix
	for (int i = dimension; i < Dim + 1; i++)
	{
		CurT[i][i] = 1;
	}

	for (int i = 0; i < dimension; i++)
	{
		CurT[i][dimension - 1] = d_axis.x[i];
		CurT[i][Dim] = points[t]->x[i];
	}

	//find those d-1 basis
	double **A, *b, *x;
	if (dimension > 1)
	{
		A = new double*[dimension - 1];
		for (int i = 0; i < dimension - 1; i++)
		{
			A[i] = new double[dimension - 1];
		}
		b = new double[dimension - 1];
		x = new double[dimension - 1];
	}
	for (int number_variables = 1; number_variables <= dimension - 1; number_variables++)
	{
		// there are number_variables variables
		//set the previous dimension-i variables to 1
		//form the equation set
		for (int i = 0; i < number_variables; i++)
		{
			double sum = 0;
			for (int j = 0; j < dimension - number_variables; j++)
			{
				if (fabs(1 - points[t]->x[j]) < ZERO_ERROR)   //when this happens, the CurT would be irreversible
				{
					sum += CurT[j][dimension - 1 - i] * (2 - points[t]->x[j]);
				}
				else
				{
					sum += CurT[j][dimension - 1 - i] * (1 - points[t]->x[j]);
				}
			}
			for (int j = dimension - number_variables; j < dimension; j++)
			{
				A[i][j - (dimension - number_variables)] = CurT[j][dimension - 1 - i];
				sum -= CurT[j][dimension - 1 - i] * points[t]->x[j];
			}
			b[i] = -sum;
		}
		bool solve = GaussianEquation(A, b, x, number_variables);
		if (!solve)
		{
			puts("the equation set is insouble!");
			break;
		}

		for (int j = 0; j < dimension - number_variables; j++)
		{
			if (fabs(1 - points[t]->x[j]) < ZERO_ERROR)   //when this happens, the CurT would be irreversible
			{
				CurT[j][dimension - 1 - number_variables] = 2 - points[t]->x[j];
			}
			else
			{
				CurT[j][dimension - 1 - number_variables] = 1 - points[t]->x[j];
			}

		}
		for (int j = dimension - number_variables; j < dimension; j++)
		{
			CurT[j][dimension - 1 - number_variables] = x[j - (dimension - number_variables)] - points[t]->x[j];
		}
	}


	//Print Current basis
	//puts("current basis ...");
	//PrintMatrix(CurT, realDimension+1);

	//compute the inverse Matrix of cur_T
	double **Inverse_M = new double*[Dim + 1];
	//store those d basis vectors
	for (int i = 0; i < Dim + 1; i++)
	{
		Inverse_M[i] = new double[Dim + 1];
	}
	bool re = GaussianInverseMatrix(CurT, Inverse_M, Dim + 1);
	if (!re)
	{
		puts("there is no inverse Matrix about CurT");   //done
		//return false;
	}
	//puts("inverse Matrix...");
	//PrintMatrix(Inverse_M, realDimension+1);


	TransformingPoints(points, Inverse_M, Dim); ////!!!!!!
	//puts("after transforming");
	//PrintPoints(points, realDimension);
	MatrixMultiply(Inverse_M, MainTainT, MainTainT, Dim + 1);
	//puts("Maintain Matrix...");
	//PrintMatrix(MainTainT, realDimension+1);

	//find the point with the largest value in current dimension
	double maxValue = 0;
	double max_index = t;
	for (int i = 0; i < n; i++)
	{
		if (points[i]->x[dimension - 1] > maxValue)
		{
			maxValue = points[i]->x[dimension - 1];
			max_index = i;
		}
	}
	curbox.L[dimension - 1] = points[t]->x[dimension - 1];
	curbox.U[dimension - 1] = points[max_index]->x[dimension - 1];

	delete[]d_axis.x;

	for (int i = 0; i < Dim + 1; i++)
	{
		delete[]CurT[i];
		delete[]Inverse_M[i];
	}
	if (dimension > 1)
	{
		for (int i = 0; i < dimension - 1; i++)
		{
			delete[]A[i];
		}
		delete[]A;
		delete[]b;
		delete[]x;
	}
	delete[]CurT;
	delete[]Inverse_M;
	return RecursionMinimumBoudingBox(points, curbox, MainTainT, dimension - 1);
}



void PrintBoudingBox(BoundingBox box)
{
	for (int i = 0; i < box.d; i++)
	{
		printf("%d dimension: %lf %lf\n", i, box.L[i], box.U[i]);
	}
}

PointSet SimpleCoreSet(PointSet &points, double **MainTainT, double epsilon)
{

	BoundingBox box(Dim);


	bool found = RecursionMinimumBoudingBox(points, box, MainTainT, Dim);
	if (!found)
	{

	}
	//PrintBoudingBox(box);
	//puts("points after computing bounding box...");
	//PrintPoints(points, dimension);

	//map the bounding box into a hypercube with length 2, centering origin
	double **HypercubeM = new double*[Dim + 1];
	for (int i = 0; i < Dim; i++)
	{
		HypercubeM[i] = new double[Dim + 1];

		for (int j = 0; j < Dim + 1; j++)
		{
			if (j < Dim)
			{
				if (j == i)
				{
					if (box.U[j] > ZERO_ERROR)
					{
						HypercubeM[i][j] = 2 / box.U[j];
					}
					else
					{
						HypercubeM[i][j] = 1;
						// non zero to make sure the matrix is reversible
					}
				}
				else
				{
					HypercubeM[i][j] = 0;
				}
			}
			else
			{
				HypercubeM[i][j] = -1;
			}

		}
	}
	HypercubeM[Dim] = new double[Dim + 1];
	for (int j = 0; j < Dim + 1; j++)
	{
		if (j < Dim)
		{
			HypercubeM[Dim][j] = 0;
			box.U[j] = 1;
			box.L[j] = -1;
		}
		else
		{
			HypercubeM[Dim][j] = 1;
		}
	}

	TransformingPoints(points, HypercubeM, Dim);
	MatrixMultiply(HypercubeM, MainTainT, MainTainT, Dim + 1);

	//puts("points after in hypercube...");
	//PrintPoints(points, dimension);
	//divided the box into pillars find the highest and lowest point in each pillar
	double C_d = 1.0 / (Dim*(4 * Dim + 1)); //a parameter
	int M = ceil(4 / (epsilon*C_d)); //number of intervals in each dimension
	if (M > 100) M = 100;
	//printf("M: %d \n", M);
	int n = points.size();
	int prime = LargestPrimeBelow(1.5*n);
	HashTable table(1.5*n, prime, Dim - 1);
	Key key(Dim - 1, M);
	for (int i = 0; i < n; i++)
	{
		//cout << " " << i << " ";
		for (int j = 0; j < Dim - 1; j++)
		{
			int cur_k = floor(1.0*(points[i]->x[j] + 1)*M / 2);
			if (cur_k >= M)
			{
				cur_k = M - 1;
			}
			if (cur_k < 0)
			{
				cur_k = 0;
			}
			key.x[j] = cur_k;
			//key = (key * M + cur_k) % prime;    //there is a problem with hash table!!!
		}
		//printf("Insert point %d\n", i);
		table.Insert(*points[i], Dim, i, key);
	}
	vector<int> index_coreset = table.Travel();

	PointSet simpleSet;
	for (int i = 0; i < index_coreset.size(); i++)
	{
		simpleSet.push_back(points[index_coreset[i]]);
	}

	for (int i = 0; i < Dim + 1; i++)
	{
		delete[]HypercubeM[i];
	}
	delete[]HypercubeM;

	return simpleSet;
}

PointSet SmallerCoreSet(PointSet &points, PointSet &directionPoints, double epsilon)
{

	double **MainTainT = new double*[Dim + 1];
	//store those d basis vectors
	for (int i = 0; i < Dim + 1; i++)
	{
		MainTainT[i] = new double[Dim + 1];
		for (int j = 0; j < Dim + 1; j++)
		{
			MainTainT[i][j] = 0;
		}
	}
	for (int i = 0; i < Dim + 1; i++)
	{
		MainTainT[i][i] = 1;
	}

	PointSet simpleSet = SimpleCoreSet(points, MainTainT, epsilon / 2);
	PointSet smallerSet;
	set<int> index_set;  //record those points added in smaller coreset
	//printf("after the simple set: size %d\n", simpleSet.size());
	//after this step, the order is disturbed
	for (int i = 0; i < directionPoints.size(); i++)
	{
		double cur_dis = MAX_DOUBLE;
		int cur_index = -1;
		for (int j = 0; j < simpleSet.size(); j++)
		{
			double tmp_dis = Distance(*directionPoints[i], *simpleSet[j], Dim);
			if (tmp_dis < cur_dis)
			{
				cur_dis = tmp_dis;
				cur_index = j;
			}
		}
		if (cur_index != -1)
		{
			if (!index_set.count(cur_index))
			{
				index_set.insert(cur_index);
				smallerSet.push_back(simpleSet[cur_index]);
			}
		}
	}

	//transforming those points into original coordinate system
	double **InverseT = new double*[Dim + 1];
	for (int i = 0; i < Dim + 1; i++)
	{
		InverseT[i] = new double[Dim + 1];
	}
	GaussianInverseMatrix(MainTainT, InverseT, Dim + 1);
	//printf("before transform\n");
	//PrintPoints(smallerSet, dimension);


	TransformingPoints(smallerSet, InverseT, Dim);

	for (int i = 0; i < Dim + 1; i++)
	{
		delete[]MainTainT[i];
		delete[]InverseT[i];
	}
	delete[]MainTainT;
	delete[]InverseT;

	return smallerSet;
}


bool DirectionalWidth(PointSet &points, HyperPlane &optimal_plane)
{
	//puts("separating positive points and negative points");
	//find +1 points and -1 points
	PointSet positivePoints, negativePoints;
	int n = points.size();
	for (int i = 0; i < n; i++)
	{
		double *x = new double[Dim];
		Point *pt = new Point(Dim, x);
		for (int j = 0; j < Dim; j++) {
			pt->x[j] = points[i]->x[j];
		}
		pt->y = points[i]->y;
		if (points[i]->y == 1)
		{
			positivePoints.push_back(pt);

		}
		else if (points[i]->y == -1)
		{
			negativePoints.push_back(pt);
		}

	}
	if (positivePoints.size() == 0 || negativePoints.size() == 0)
	{
		puts("thre are only one class in the point set");
		return false;
	}
	//compute those directions set the angle with epsilon
	//compute direction set for constructing CoreSet
	double alpha = 1.0 / (Dim*(4 * Dim + 1));
	double delta = sqrt(Rho*alpha / 4);
	double radius = sqrt(Dim) + 1;


	double *angles = new double[Dim - 1];
	PointSet directionPoints;
	//puts("computing directions...");
	ComputingDirections(directionPoints, angles, 1, delta, radius, 2*M_PI);

	PointSet positiveCoreSet = SmallerCoreSet(positivePoints, directionPoints, Rho);
	PointSet negativeCoreSet = SmallerCoreSet(negativePoints, directionPoints, Rho);

	//merge the two core set
	positiveCoreSet.insert(positiveCoreSet.end(), negativeCoreSet.begin(), negativeCoreSet.end());
	
	//compute a hyperplane along each direction
	//actually the angle of core set and classification may be different !!!!!
	//a different direction points the angle is epsilon
	radius = 1;
	directionPoints.clear();
	ComputingDirections(directionPoints, angles, 1, Rho, radius, M_PI);
	//printf("there are %d direction points\n", directionPoints.size());
	double max_margin = 0;
	//
	for (int i = 0; i < directionPoints.size(); i++)
	{
		HyperPlane cur_plane(Dim);
		bool found = OneDimensionClassification(positiveCoreSet, cur_plane, *directionPoints[i], radius);
		//compute the margin of current plane
		if (found)
		{
			//puts("separable on those coreset points");
			//PrintHyperPlane(cur_plane, dimension);
			//check each hyperplane direction with the original data
			double cur_distance = MAX_DOUBLE;
			bool separable = MinimumSeparableDistance(points, cur_plane, cur_distance);
			if (separable && cur_distance > max_margin)
			{
				//optimal_plane = cur_plane; // cannot assign directly
				CopyHyperPlane(optimal_plane, cur_plane);
				max_margin = cur_distance;
			}
		}
	}

	delete[]angles;
	int posSize = positivePoints.size();
	int negSize = negativePoints.size();
	for (int i = 0; i < posSize; i++) {
		delete []positivePoints[i]->x;
		delete (positivePoints[i]);
		positivePoints[i] = NULL;
	}
	for (int i = 0; i < negSize; i++) {
		delete[]negativePoints[i]->x;
		delete (negativePoints[i]);
		negativePoints[i] = NULL;
	}

	//return the largest margin hyperplane
	if (max_margin > ZERO)
	{
		PrintHyperPlane(optimal_plane);
		cout << max_margin << endl;
		puts("find an solution along those directions");
		return true;
	}
	else
	{
		return false;
	}

}

bool DirectionalWidth(PointSet &points, HyperPlane &optimal_plane,
		PointSet &coresetDirections, PointSet &classifyDirections)
{
	//puts("separating positive points and negative points");
	//find +1 points and -1 points
	PointSet positivePoints, negativePoints;
	int n = points.size();
	for (int i = 0; i < n; i++)
	{
		double *x = new double[Dim];
		Point *pt = new Point(Dim, x);
		for (int j = 0; j < Dim; j++) {
			pt->x[j] = points[i]->x[j];
		}
		pt->y = points[i]->y;
		if (points[i]->y == 1)
		{
			positivePoints.push_back(pt);

		}
		else if (points[i]->y == -1)
		{
			negativePoints.push_back(pt);
		}

	}
	if (positivePoints.size() == 0 || negativePoints.size() == 0)
	{
		puts("thre are only one class in the point set");
		return false;
	}
	puts("computing coreset");
	//cout << positivePoints.size() << " " << negativePoints.size() << endl;
	PointSet positiveCoreSet = SmallerCoreSet(positivePoints, coresetDirections, Rho);
	PointSet negativeCoreSet = SmallerCoreSet(negativePoints, coresetDirections, Rho);
	
	//merge the two core set
	positiveCoreSet.insert(positiveCoreSet.end(), negativeCoreSet.begin(), negativeCoreSet.end());
	//printf("after the smaller set: size %d\n", positiveCoreSet.size());
	puts("computing classifier");
	double max_margin = 0;
	//
	for (int i = 0; i < classifyDirections.size(); i++)
	{
		HyperPlane cur_plane(Dim);
		bool found = OneDimensionClassification(positiveCoreSet, cur_plane, *classifyDirections[i], 1);
		//compute the margin of current plane
		if (found)
		{
			//puts("separable on those coreset points");
			//PrintHyperPlane(cur_plane, dimension);
			//check each hyperplane direction with the original data
			double cur_distance = MAX_DOUBLE;
			bool separable = MinimumSeparableDistance(points, cur_plane, cur_distance);
			if (separable && cur_distance > max_margin)
			{
				//optimal_plane = cur_plane; // cannot assign directly
				CopyHyperPlane(optimal_plane, cur_plane);
				max_margin = cur_distance;
			}
		}
	}
	int posSize = positivePoints.size();
	int negSize = negativePoints.size();
	for (int i = 0; i < posSize; i++) {
		delete[]positivePoints[i]->x;
		delete (positivePoints[i]);
		positivePoints[i] = NULL;
	}
	for (int i = 0; i < negSize; i++) {
		delete[]negativePoints[i]->x;
		delete (negativePoints[i]);
		negativePoints[i] = NULL;
	}
	//return the largest margin hyperplane
	if (max_margin > ZERO)
	{
		PrintHyperPlane(optimal_plane);
		cout << max_margin << endl;
		puts("find an solution along those directions");
		return true;
	}
	else
	{
		return false;
	}
}

bool OneDimensionClassification(PointSet &points, HyperPlane &cur_plane, Point &direction, double radius)
{
	int n = points.size();
	OneDim *one_pts = new OneDim[n];
	for (int i = 0; i < n; i++)
	{
		one_pts[i].x = Dot(*points[i], direction, cur_plane.d) / radius;
		one_pts[i].y = points[i]->y;
	}
	sort(one_pts, one_pts + n);
	int *prefix_sum_positive = new int[n];
	int *suffix_sum_negative = new int[n];

	for (int i = 0; i < n; i++)
		suffix_sum_negative[i] = 0;
	for (int i = 0; i < n; i++)
		prefix_sum_positive[i] = 0;
	if (one_pts[0].y == 1)
	{
		prefix_sum_positive[0] = 1;
	}
	else
	{
		prefix_sum_positive[0] = 0;
	}
	for (int i = 1; i < n; i++)
	{
		if (one_pts[i].y == 1)
			prefix_sum_positive[i] = prefix_sum_positive[i - 1] + 1;
		else
			prefix_sum_positive[i] = prefix_sum_positive[i - 1];
	}

	if (one_pts[n - 1].y == -1)
	{
		suffix_sum_negative[n - 1] = 1;
	}
	else
	{
		suffix_sum_negative[n - 1] = 0;
	}
	for (int i = n - 2; i >= 0; i--)
	{
		if (one_pts[i].y == -1)
			suffix_sum_negative[i] = suffix_sum_negative[i + 1] + 1;
		else
			suffix_sum_negative[i] = suffix_sum_negative[i + 1];
	}
	bool separable = false;
	double line = 0;
	for (int i = 0; i < n; i++)
	{
		int total = prefix_sum_positive[i] + suffix_sum_negative[i];
		if (total == n || total == 1)
		{
			separable = true;
			line = (one_pts[i].x + one_pts[i + 1].x) / 2;
			break;

		}
	}
	delete[]prefix_sum_positive;
	delete[]suffix_sum_negative;
	delete[]one_pts;

	if (!separable)
	{
		return false;
	}
	int dimension = cur_plane.d;
	double sum = 0;
	for (int j = 0; j < dimension; j++)
	{
		cur_plane.w[j] = direction.x[j] / radius;
		sum = sum - line*(direction.x[j] / radius) * cur_plane.w[j];
	}
	cur_plane.b = sum;

	return separable;

}

void ComputingDirections(PointSet &points, double*angles, int curDimension, double delta, double radius, double full_angle)
{
	if (curDimension == Dim)
	{
		double *x = new double[Dim];
		Point *cur_pt = new Point(Dim, x);
		for (int j = 0; j < Dim; j++)
		{
			cur_pt->x[j] = radius;
		}
		for (int i = 0; i < Dim - 1; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				if (j == i)
				{
					cur_pt->x[i] *= sin(angles[j]);
				}
				else
				{
					cur_pt->x[i] *= cos(angles[j]);
				}
			}
		}
		for (int j = 0; j < Dim - 1; j++)
		{
			cur_pt->x[Dim - 1] *= cos(angles[j]);
		}
		points.push_back(cur_pt);
		/*
		if (points.size() % 10000 == 0)
		{
		cout << "current number of direcitons: ";
		cout << points.size() << endl;
		}
		*/
		return;
	}
	double cur_angle = 0;
	for (; cur_angle < full_angle; cur_angle += delta)
	{
		angles[curDimension - 1] = cur_angle;
		ComputingDirections(points, angles, curDimension + 1, delta, radius, full_angle);
	}
}

