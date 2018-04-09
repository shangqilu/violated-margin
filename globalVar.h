#ifndef GLOBALVAR_H
#define GLOBALVAR_H

/*
*  the number of points
*/
extern int N;


/*
*  the dimensionality of points
*/
extern int Dim;

/*
*  the largest classification error rate accepted
*/
extern double Percent;


/*
*  the largest number of points classifid wrongly accepted
*/
extern int K;

/*
*  the approximation for error rate
*/
extern double Epsilon;

/*
*  the approximation for margin
*/
extern double Rho;
/*
*  The failing probability
*/
extern double Delta;
/*
*  The chosen method.
*/
extern int Method;





#endif //GLOBALVAR_H