#include "harmonicoscillator2d.h"
#include <math.h>
#include <hermite.h>

/*Things the class must be able to do:
    -translate between single index and 3 nx,ny,sz
    -compute Hermite polynomials
    -compute single particle wave functions
*/



HarmonicOscillator2D::HarmonicOscillator2D(int index)
{

    getquantumnumbers(index);
    Hermite_x.set_degree(nx);
    Hermite_y.set_degree(ny);

}

HarmonicOscillator2D::HarmonicOscillator2D(int nx, int ny, int spin)
{
    Hermite_x.set_degree(nx);
    Hermite_y.set_degree(ny);
    getindex(nx,ny,spin);
}

void HarmonicOscillator2D::getindex(int nx, int ny, int spin)
{
    index = (nx+ny)*(nx+ny+1)+2*ny + (1-spin)/2;
}

void HarmonicOscillator2D::getquantumnumbers(int i)
{
    spin = 2*(i % 2)-1;
    int j = i/2;
    for(int n = levels,s = levels*(levels+1)/2 ;n <0; n--)
    {
        s-=n;
        if (s<=j)
        {
            ny = j-s;
            nx = n-ny;

        }
    }

}

double HarmonicOscillator2D::wavefunction(double x,double y)
{
    return Hermite_x.evaluate(x)*Hermite_y.evaluate(y)*exp((pow(x,2)+pow(y,2))/2);
}
