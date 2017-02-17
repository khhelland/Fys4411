#ifndef HARMONICOSCILLATOR2D_H
#define HARMONICOSCILLATOR2D_H
#include "hermite.h"

class HarmonicOscillator2D
{
public:
    HarmonicOscillator2D();
    HarmonicOscillator2D(int);
    HarmonicOscillator2D(int index, double w);
//    HarmonicOscillator2D(int, int, int, double w);

//    double wavefunction_cartesian(double, double y);
//    double wavefunction_cartesian_no_exp(double,double);

//    int getn();
//    int getm();

    int n;
    int m;

private:
    int levels = 10;
    int level;
    int spin;
    int index;

    double omega = 1;

//    int nx;
//    int ny;
//    double normalization_constant;

//    double pi = 3.141592653589793238462643383279502884;

//    Hermite Hermite_x;
//    Hermite Hermite_y;

//    void setindex(int nx, int ny, int spin);
//    void setquantumnumbers_cartesian(int i);
    void setquantumnumbers_polar(int i);
    void setlevel(int);


//    int factorial(int x);
//    void normalize();




};

#endif // HARMONICOSCILLATOR2D_H
