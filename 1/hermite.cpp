#include "hermite.h"
#include <math.h>
#include <iostream>

using namespace std;

Hermite::Hermite(int deg)
{
    degree = deg;

}

Hermite::Hermite(){}

void Hermite::set_degree(int deg)
{
    degree = deg;
}

double Hermite::evaluate(double x)
{
    double result;
    switch(degree)
    {
        case 0 : result = 1;
                 break;
        case 1 : result = 2*x;
                 break;
        case 2 : result = 4*pow(x,2)-2;
                 break;
        case 3 : result = 8*pow(x,3)-12*x;
                 break;
        case 4 : result = 16*pow(x,4)-48*pow(x,2) + 12;
                 break;
        case 5 : result = 32*pow(x,5) - 160*pow(x,3) + 120*x;
                 break;
        case 6 : result = 64*pow(x,6) - 480*pow(x,4) + 720*pow(x,2)-120;
                 break;
        case 7 : result = 128*pow(x,7) - 1344*pow(x,5) + 3360*pow(x,3) - 1680*x;
                 break;
        case 8 : result = 256*pow(x,8) - 3584*pow(x,6) + 13440*pow(x,4) - 13440*pow(x,2) + 1680;
                 break;
        case 9 : result = 512*pow(x,9) -9216*pow(x,7) + 48384*pow(x,5) -80640*pow(x,3) + 30240*x;
    }
    return result;
}
