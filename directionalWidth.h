#ifndef DIRECTIONALWIDTH_H_INCLUDED
#define DIRECTIONALWIDTH_H_INCLUDED

#include "headers.h"

using namespace std;

//max iterations to find the two approximation diameter
const int MaxIterTwoApprox = 1000;

/*
*   structure of the Bounding Box
*/
class BoundingBox{
public:
	double *U;
	double *L;
	int d;

	BoundingBox(int d)
	{
		this->U = new double[d];
		this->L = new double[d];
		this->d = d;
	}
	~BoundingBox()
	{
		if (U != NULL) {
			delete[]U;
			U = NULL;
		}
		if (L != NULL) {
			delete[]L;
			L = NULL;
		}
	}
private:
	BoundingBox(const BoundingBox &obj){}
	BoundingBox & operator = (const BoundingBox & obj){}
};

/*
*   structure of a key in hash table
*/
class Key
{
public:
	int *x;
	int d;
	int base;
	Key(int d, int base)
	{
		this->d = d;
		this->base = base;
		this->x = new int[d];
	}
	~Key()
	{
		if (x != NULL) {
			delete[]x;
			x = NULL;
		}
	}
private:
	Key(const Key &obj){}
	Key & operator = (const Key &obj){}
};

/*
*   structure of a pair in hash table
*/
class Pair{
public:
	int low_index;      //the index of lowest point in current pillar
	int high_index;     //the index of highest point in current pillar
	double low_value;   //the value of lowest point in current pillar
	double high_value;  //the value of lowest point in current pillar
	int *x;             //the coordinate of current pillar
	int d;
	Pair(int d)
	{
		this->x = new int[d];
		this->d = d;
		low_index = high_index = -1;
		low_value = high_value = 0;
	}
	~Pair()
	{
		if (x != NULL) {
			delete[]x;
			x = NULL;
		}
	}
	Pair(const Pair &obj)
	{
		x = new int[obj.d];
		for (int i = 0; i < obj.d; i++)
		{
			x[i] = obj.x[i];
		}
		low_index = obj.low_index;
		high_index = obj.high_index;
		low_value = obj.low_value;
		high_value = obj.high_value;
		d = obj.d;
	}
	Pair& operator = (const Pair &obj)
	{
		if (this == &obj) {
			return *this;
		}
		else {
			delete[]x;
			x = new int[obj.d];
			for (int i = 0; i < obj.d; i++) {
				x[i] = obj.x[i];
			}
			low_index = obj.low_index;
			high_index = obj.high_index;
			low_value = obj.low_value;
			high_value = obj.high_value;
			d = obj.d;
		}
	}
private:
	
};

/*
*   structure to represent a point in one dimension
*/
struct OneDim{
    double x;
    int y;
    bool operator < (const OneDim & a) const
    {
        return this->x < a.x;
    }
};

/*
*   hash table to find the lowest and highest points in a pillar
*/
class HashTable{
public:
	vector<Pair> table;
	bool *Empty;    //flag array
	int M;          //table size
	int p;          //divisor
	HashTable(int M, int p, int d) //using the division method
	{
		this->M = M;
		this->p = p;
		for (int i = 0; i < M; i++) {
			Pair newpair(d);
			table.push_back(newpair);
		}
		Empty = new bool[M];
		for (int i = 0; i < M; i++) {
			Empty[i] = true;
		}
	}
	~HashTable()
	{
		if (Empty != NULL) {
			delete[]Empty;
			Empty = NULL;
		}
	}

    int Hash_func(Key &key)
    {
        int value = 0;
        for (int i = 0; i < key.d; i++) {
            value = ((value%p) * (key.base%p) + key.x[i]) % p;
        }
        return value;
    }

    bool key_match(Pair cur_pair, Key &key)
    {
        for (int i = 0; i < key.d; i++) {
            if (cur_pair.x[i] != key.x[i]) return false;
        }
        return true;
    }
    void PrintKey(Key &key)
    {
        cout << "(";
        for (int i = 0; i < key.d-1; i++) {
            cout << key.x[i] << ", ";
        }
        cout << key.x[key.d-1];
        cout << ")" << endl;
    }
    bool Insert(Point &point, int d, int index, Key &key)
    {
        int value = Hash_func(key);
        value = value - 1;
        int di = (value + 1) % M;
        while(1) {
            if (Empty[di]) {
                break;
            }
            if (!Empty[di] && key_match(table[di], key)) {
                break;
            }
            if (value == -1 && (di+1) == M) {
                return false;
            }
            if (di != value) {
                di = (di + 1) % M;
            }else {
                return false;
            }
        }
        //insert current point
        if (Empty[di]) {
            //printf("Empty: index %d ", index);
            //PrintKey(key);
            Empty[di] = false;
            table[di].high_index = index;
            table[di].high_value = point.x[d-1];
            table[di].low_index = index;
            table[di].low_value = point.x[d-1];
            for (int i = 0; i < key.d; i++) {
                table[di].x[i] = key.x[i];
            }
            return true;
        }
        //find the slot and update
        if (!Empty[di] && key_match(table[di], key)) {
            //printf("Update: index %d ", index);
            //PrintKey(key);
            if (point.x[d-1] > table[di].high_value) {
                table[di].high_value = point.x[d-1];
                table[di].high_index = index;
            }
            if (point.x[d-1] < table[di].low_value) {
                table[di].low_value = point.x[d-1];
                table[di].low_index = index;
            }
            return true;
        }
		return false;
    }

    vector<int> Travel() //report all the points in hash table
    {
        vector<int> index_points;
        //puts("traveling the hash table");
        for (int i = 0; i < M; i++)
        {
            if (!Empty[i]) {
                index_points.push_back(table[i].high_index);
                if (table[i].low_index != table[i].high_index)
                {
                    index_points.push_back(table[i].low_index);
                }
            }
        }
        //cout << index_points.size() << endl;
        return index_points;
    }
private:
	HashTable& operator = (const HashTable &obj){}
	HashTable(const HashTable &obj){}
};

void PrintBoundingBox(BoundingBox box);


/*
*   Compute the bounding box recursively
*   Parameter List:
*       points: points set
*       curbox: current bounding box
*       MainTainT: the accumulated transforming
*       dimension: current dimension
*       realDimension: the real dimension of points and bounding box
*/
bool RecursionMinimumBoudingBox(PointSet &points, BoundingBox &curbox, double **MainTainT, int dimension, int realDimension);

/*
*   Compute the two approximate diameter of a point set
*   Parameter List:
*       points: points set
*       s: the index of a point
*       t: the index of the other point, the other point lies one side of the plane
*           whose normal vector is ts and point t is on the plane
*/
bool TwoApproxiDiameter(PointSet &points, int dimension, int &s, int &t);

/*
*   compute a simple core set with size O(1/epsilon^d-1)
*/
PointSet SimpleCoreSet(PointSet &points, double **MainTainT, int dimension, double epsilon);

/*
*   compute a smaller core set with size O(1/epsilon^(d-1)/2)
*/
PointSet SmallerCoreSet(PointSet &points, PointSet directionPoints, int dimension, double epsilon);

/*
*   Compute a hyperplane with a (1-rho)approximation margin using Directional width algorithm
*/
bool DirectionalWidth(PointSet &points, HyperPlane &plane, int dimension, double rho);

/*
*   Given a direction, compute a hyperplane with largest margin in one dimension
*/
bool OneDimensionClassification(PointSet &points, HyperPlane &cur_plane, Point direction, double radius);

/*
*   recursively compute directions
*   store them in points set
*/
void ComputingDirections(PointSet &points, double*angles, int curDimension, int realDimension, double delta, double radius);


#endif // DIRECTIONALWIDTH_H_INCLUDED
