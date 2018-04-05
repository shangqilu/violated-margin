#include "test.h"
#include "headers.h"
#include <iostream>

using namespace std;




int main(int nargs, char **args)
{
    /*
    if (nargs != 2 || !isdigit(args[1][0])) {
        puts("input: executing_name method_number");
        puts("there are three methods");
        puts("0: perceptron");
        puts("1: simplex");
        puts("2: directional width");
        return 0;
    }
    */
    //int method = args[1][0] - '0';
    //cout << method << endl;
    //TestViolatedMargin(method);
    
    //TestDirectionalWidth();
    //TestPerceptron();
    TestSimplex();
    return 0;
}
