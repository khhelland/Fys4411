#ifndef HARMONICOSCILLATOR2D_H
#define HARMONICOSCILLATOR2D_H
#include "hermite.h"

class HarmonicOscillator2D
{
public:
    HarmonicOscillator2D(int index, double w);
    HarmonicOscillator2D(int, int, int, double w);

    double wavefunction(double, double y);
    double wavefunction_no_exp(double,double);

private:
    int levels = 10;
    int spin;
    int nx;
    int ny;
    int index;
    double normalization_constant;
    double omega;
    const double pi = 3.141592653589793238462643383279502884;

    Hermite Hermite_x;
    Hermite Hermite_y;

    void getindex(int nx, int ny, int spin);
    void getquantumnumbers(int i);
    int factorial(int x);
    void normalize();




};

#endif // HARMONICOSCILLATOR2D_H
