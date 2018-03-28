#ifndef DIRECTIONALWIDTH_H_INCLUDED
#define DIRECTIONALWIDTH_H_INCLUDED

#include "headers.h"

using namespace std;

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

struct Pair{
    int low_index; //the index of lowest point in current pillar
    int high_index; //the index of highest point in current pillar
    double low_value;
    double high_value;
    int key;
};

struct HashTable{
    Pair *table;
    bool *Empty;
    int M; //table size
    int p; //divisor
    HashTable(int M, int p) //using the division method
    {
        this->M = M;
        this->p = p;
        table = new Pair[M];
        Empty = new bool[M];
        for (int i = 0; i < M; i++) {
            Empty[i] = true;
        }
    }

    int Hash_func(int key)
    {
        return key % p;
    }


    bool Insert(Point point, int d, int index, int key)
    {
        int value = Hash_func(key);
        value = value - 1;
        int di = (value + 1) % M;
        while(!Empty[di] && table[di].key != key) {
            if (di != value) {
                di = (di + 1) % M;
            }else {
                return false;
            }
        }
        if (Empty[di]) {  //insert current point
            Empty[di] = false;
            table[di].high_index = index;
            table[di].high_value = point.x[d-1];
            table[di].low_index = index;
            table[di].low_value = point.x[d-1];
            return true;
        }

        if (!Empty[di] && key == table[di].key) {//find the slot and update
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
        for (int i = 0; i < M; i++)
        {
            if (!Empty[i]) {
                index_points.push_back(table[i].high_index);
                index_points.push_back(table[i].low_index);
            }
        }
        return index_points;
    }

};

void PrintBoudingBox(BoudingBox box);
void MinimumBoudingBox(PointSet points, int dimension);

void RecursionMinimumBoudingBox(PointSet points, BoudingBox &curbox, double **MainTainT, int dimension, int realDimension);
void TwoApproxiDiameter(PointSet points, int dimension, int &s, int &t);

#endif // DIRECTIONALWIDTH_H_INCLUDED
