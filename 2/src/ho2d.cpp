#include "ho2d.h"
#include <iostream>
#include <math.h>

using namespace std;

/*
 * Using case instead of finding formulas.
 * This limits the usage to 20 particles.
 * If extending try generalizing.
 */

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

double ho2dDiff(int deg,double w, double x, int dim)
{
    int nx = 0;
    int ny = 0;
    switch (deg) {
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
    int n = (dim==0)? nx : ny;
    return (sq*hermitederiv(n,sq*x)-x*hermite(n,sq*x))*exp(-x*x/2);
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

double hermitederiv(int d, double x)
{
    double result = 0;
    switch(d)
    {
    case 0 : result = 0;
        break;
    case 1 : result = 2;
        break;
    case 2 : result = 8*x;
        break;
    case 3 : result = 24*x*x-12;
        break;
    case 4 : result = 64*x*x*x-96*x ;
        break;
    default:
        cout<<"Warning: derivative of Hermite polynomial not implemented"<<endl;
    }
    return result;
}



int ho2denergy(int nOrbitals)
{
    //energy of nOrbitals lowest orbitals
    int result = 0;
    switch(nOrbitals)
    {
    case 1 : result = 1;
        break;
    case 3 : result = 5;
        break;
    case 6 : result = 10;
        break;
    case 10 : result = 30;
        break;
    default:
        cout<<"Warning: ground state energy not implemented for this N"<<endl;
    }
    return result;
}


