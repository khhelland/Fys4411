#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <harmonicoscillator2d.h>


class integrator
{
public:
    integrator(HarmonicOscillator2D,HarmonicOscillator2D,HarmonicOscillator2D,HarmonicOscillator2D);
    double integrate(int n);

private:
    HarmonicOscillator2D p;
    HarmonicOscillator2D q;
    HarmonicOscillator2D r;
    HarmonicOscillator2D s;

    // The integrand in Cartesian coordinates
    double Integrand(double, double,double, double);
    // The Gauss Hermite integration function
    double  GaussHermiteIntegration(int);

    // Getting the Gaussian quadrature weights and integration points
    void GaussHermiteQuadrature(double *, double *, int);


};

#endif // INTEGRATOR_H
