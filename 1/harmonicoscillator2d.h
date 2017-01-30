#ifndef HARMONICOSCILLATOR2D_H
#define HARMONICOSCILLATOR2D_H
#include "hermite.h"

class HarmonicOscillator2D
{
public:
    HarmonicOscillator2D(int index);
    HarmonicOscillator2D(int,int,int);

    double wavefunction(double, double y);

private:
    int levels = 10;
    int spin;
    int nx;
    int ny;
    int index;
    Hermite Hermite_x;
    Hermite Hermite_y;

    void getindex(int nx, int ny, int spin);
    void getquantumnumbers(int i);





};

#endif // HARMONICOSCILLATOR2D_H
