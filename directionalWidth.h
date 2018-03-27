#ifndef DIRECTIONALWIDTH_H_INCLUDED
#define DIRECTIONALWIDTH_H_INCLUDED

#include "headers.h"


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




#endif // DIRECTIONALWIDTH_H_INCLUDED
