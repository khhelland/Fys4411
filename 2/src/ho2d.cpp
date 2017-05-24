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
    arma::vec ns = quantumnumbers(d);
    int nx = ns(0);
    int ny = ns(1);
    double sq = sqrt(w);

    return hermite(nx,sq*x)*hermite(ny,sq*y)*exp(-w*(x*x + y*y)/2);
}

arma::vec ho2dgrad(int deg, double w, double x, double y)
{
    arma::vec ns = quantumnumbers(deg);
    int nx = ns(0);
    int ny = ns(1);

    double sq = sqrt(w);

    arma::vec grad(2);

    grad(0) = (sq*hermitederiv(nx,sq*x,1)-x*w*hermite(nx,sq*x))*hermite(ny,sq*y);
    grad(1) = (sq*hermitederiv(ny,sq*y,1)-y*w*hermite(ny,sq*y))*hermite(nx,sq*x);
    grad *= exp(-w*(x*x + y*y)/2);
    return grad;

}

arma::vec quantumnumbers(int deg)
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
    arma::vec result = {nx,ny};
    return result;
}


double hermite(int d, double x)
{
    double result = 0;
    switch(d)
    {
    case -1 : result = 0;
        break;
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

double hermitederiv(int d, double x,int n)
{
    //nth derivative of hermite poly of deg d at x
    if(n>d) return 0;
    double fac = 1;
    for(int i = n; i > 0; i--) fac*=2*(d-i+1);
    return fac*hermite(d-n,x);
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
    case 6 : result = 14;
        break;
    case 10 : result = 30;
        break;
    default:
        cout<<"Warning: ground state energy not implemented for this N"<<endl;
    }
    return result;
}


