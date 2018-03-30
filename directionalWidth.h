#ifndef DIRECTIONALWIDTH_H_INCLUDED
#define DIRECTIONALWIDTH_H_INCLUDED

#include "headers.h"

using namespace std;

const int MaxIterTwoApprox = 1000;

struct BoudingBox{
    double *U;
    double *L;
    int d;

    BoudingBox(int d)
    {
        this->U = new double[d];
        this->L = new double[d];
        this->d = d;
    }

};

struct Key
{
    int *x;
    int d;
    int base;
    Key(int d, int base)
    {
        this->d = d;
        this->base = base;
        this->x = new int[d];
    }
};

struct Pair{
    int low_index; //the index of lowest point in current pillar
    int high_index; //the index of highest point in current pillar
    double low_value;
    double high_value;
    int *x;
    Pair(int d)
    {
        this->x = new int[d];
    }
};

struct OneDim{
    double x;
    int y;

    bool operator < (const OneDim & a) const
    {
        return this->x < a.x;
    }

};


struct HashTable{
    vector<Pair> table;
    bool *Empty;
    int M; //table size
    int p; //divisor
    HashTable(int M, int p, int d) //using the division method
    {
        this->M = M;
        this->p = p;
        for (int i = 0; i < M; i++) {
            Pair newpair = Pair(d);
            table.push_back(newpair);
        }
        Empty = new bool[M];
        for (int i = 0; i < M; i++) {
            Empty[i] = true;
        }
    }

    int Hash_func(Key key)
    {
        int value = 0;
        for (int i = 0; i < key.d; i++) {
            value = ((value%p) * (key.base%p) + key.x[i]) % p;
        }
        return value;
    }

    bool key_match(Pair cur_pair, Key key2)
    {
        for (int i = 0; i < key2.d; i++) {
            if (cur_pair.x[i] != key2.x[i]) return false;
        }
        return true;
    }
    void PrintKey(Key key)
    {
        cout << "(";
        for (int i = 0; i < key.d-1; i++) {
            cout << key.x[i] << ", ";
        }
        cout << key.x[key.d-1];
        cout << ")" << endl;
    }
    bool Insert(Point point, int d, int index, Key key)
    {
        int value = Hash_func(key);
        value = value - 1;
        int di = (value + 1) % M;
        while(!Empty[di] && !key_match(table[di], key)) {
            if (di != value) {
                di = (di + 1) % M;
            }else {
                return false;
            }
        }
        if (Empty[di]) {  //insert current point
            printf("Empty: index %d ", index);
            PrintKey(key);
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

        if (!Empty[di] && key_match(table[di], key)) {//find the slot and update
            printf("update: index %d", index);
            PrintKey(key);
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
    }

    vector<int> Travel() //report all the points in hash table
    {
        vector<int> index_points;
        puts("traveling the hash table");
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
        cout << index_points.size() << endl;
        return index_points;
    }

};

void PrintBoudingBox(BoudingBox box);
void MinimumBoudingBox(PointSet points, int dimension);

bool RecursionMinimumBoudingBox(PointSet &points, BoudingBox &curbox, double **MainTainT, int dimension, int realDimension);
bool TwoApproxiDiameter(PointSet points, int dimension, int &s, int &t);
PointSet SimpleCoreSet(PointSet points, double **MainTainT, int dimension, double epsilon);

PointSet SmallerCoreSet(PointSet points, PointSet directionPoints, int dimension, double epsilon);

bool DirectionalWidth(PointSet points, HyperPlane &plane, int dimension, double epsilon);

bool OneDimensionClassification(PointSet points, HyperPlane &cur_plane, Point direction, double radius);

void ComputingDirections(PointSet &points, double*angles, int curDimension, int realDimension, double delta, double radius);


#endif // DIRECTIONALWIDTH_H_INCLUDED
