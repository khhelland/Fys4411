#include "ho2d.h"
#include <iostream>
#include <math.h>

using namespace std;

double ho2d(int d, double w, double x, double y)
{
    int nx = 0;
    int ny = 0;
    switch (d) {
    case 0:
        nx = 0;
        ny = 0;
        break;
    case 1:
        nx = 1;
        ny = 0;
        break;
    case 2:
        nx = 0;
        ny = 1;
        break;
    case 3:
        nx = 2;
        ny = 0;
        break;
    case 4:
        nx = 1;
        ny = 1;
        break;
    case 5:
        nx = 0;
        ny = 2;
        break;
    case 6:
        nx = 3;
        ny = 0;
        break;
    case 7:
        nx = 2;
        ny = 1;
        break;
    case 8:
        nx = 1;
        ny = 2;
        break;
    case 9:
        nx = 0;
        ny = 3;
        break;
    default:
        cout<<"Warning: single particle wavefunction not implemented"<<endl;
        }
    double sq = sqrt(w);

    return hermite(nx,sq*x)*hermite(ny,sq*y)*exp(-w*(x*x + y*y)/2);
}


double hermite(int d, double x)
{
    double result = 0;
    switch(d)
    {
    case 0 : result = 1;
        break;
    case 1 : result = 2*x;
        break;
    case 2 : result = 4*x*x-2;
        break;
    case 3 : result = 8*x*x*x-12*x;
        break;
    case 4 : result = 16*x*x*x*x-48*x*x + 12;
        break;
    default:
        cout<<"Warning: Hermite polynomial not implemented"<<endl;
    }
    return result;
}




