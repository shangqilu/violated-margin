#ifndef TEST_H_INCLUDED
#define TEST_H_INCLUDED

#include "headers.h"
#include "perceptron.h"
#include "simplex.h"
#include "directionalWidth.h"

using namespace std;


void TestSimplex();

void TestSimplex2();

void TestPerceptron();

void TestGaussianEquation();
void TestGaussianInverse();

void TestBoudingBox();
void TestSimpleCoreSet();


void TestSmallerCoreSetandComputingDirection();

void TestOneDimensionClassification();

void TestDirectionalWidth();

#endif // TEST_H_INCLUDED
