#include "harmonicoscillator2d.h"
#include <math.h>
#include <hermite.h>
#include <iostream>
/*Things the class must be able to do:
    -translate between single index and 3 nx,ny,sz
    -compute Hermite polynomials
    -compute single particle wave functions
*/
// --------------- Constructors-------------------------------
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
HarmonicOscillator2D::HarmonicOscillator2D()
{

}

HarmonicOscillator2D::HarmonicOscillator2D(int index1)
{
    index = index1;

    setquantumnumbers_polar(index);
//    std::cout<<level<<" "<<n<<" "<<m<<std::endl;

//    setquantumnumbers_cartesian(index);
//    Hermite_x.set_degree(nx);
//    Hermite_y.set_degree(ny);
//    normalize();
}

HarmonicOscillator2D::HarmonicOscillator2D(int index1,double w)
{
    omega = w;
    index = index1;

    setquantumnumbers_polar(index);

//    setquantumnumbers(index);
//    Hermite_x.set_degree(nx);
//    Hermite_y.set_degree(ny);
//    normalize();
}

double HarmonicOscillator2D::getenergy()
{
    return omega*(2*n+abs(m)+1);
}

//HarmonicOscillator2D::HarmonicOscillator2D(int nx, int ny, int spin, double w)
//{
//    omega = w;
//    Hermite_x.set_degree(nx);
//    Hermite_y.set_degree(ny);
//    getindex(nx,ny,spin);
//}
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



//-----------------------Quantum Numbers------------------------------------------
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//void HarmonicOscillator2D::setindex(int nx, int ny, int spin)
//{
//    index = (nx+ny)*(nx+ny+1)+2*ny + (1-spin)/2;
//}

//void HarmonicOscillator2D::setquantumnumbers_cartesian(int i)
//{
//    spin = 2*(i % 2)-1;
//    int orbital = i/2;
//    setlevel(i);
//    ny = j-level(level-1)/2;
//    nx = level-1-ny;
//}

void HarmonicOscillator2D::setquantumnumbers_polar(int i)
{
    spin = 2*(i%2)-1;
    int orbital = i/2;
    setlevel(i);
    m = 2*orbital - level*level + 1;
    n = (level -abs(m) -1)/2;
}

void HarmonicOscillator2D::setlevel(int p)
{
    int orbital = p/2;
    int current_level = levels;
    for(int s = levels*(levels+1)/2; s>orbital ; current_level--)
    {
        s-=current_level;
    }
    level = current_level+1;
}


//int HarmonicOscillator2D::getm() {return m}
//int HarmonicOscillator2D::getn() {return n}

//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


////---------------------Wavefunction------------------------------------
////||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//double HarmonicOscillator2D::wavefunction_cartesian(double x,double y)
//{
//    return wavefunction_cartesian_no_exp(x,y)*exp((pow(x,2)+pow(y,2))/2);
//}

//double HarmonicOscillator2D::wavefunction_cartesian_no_exp(double x, double y)
//{
//    return normalization_constant*Hermite_x.evaluate(x)*Hermite_y.evaluate(y);
//}


//void HarmonicOscillator2D::normalize()
//{
//    normalization_constant = pow(omega/pi,0.5)/sqrt(pow(2,nx+ny)*factorial(nx)*factorial(ny));
//}

//int HarmonicOscillator2D::factorial(int x)
//{
//    return (x == 1 || x==0) ? 1 : factorial(x-1)*x;}
////|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
