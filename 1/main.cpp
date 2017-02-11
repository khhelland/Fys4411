#include <iostream>
#include "harmonicoscillator2d.h"
#include <hermite.h>
#include <integrator.h>


using namespace std;

int main()
{
    HarmonicOscillator2D p(5);
    HarmonicOscillator2D q(1);
    HarmonicOscillator2D r(4);
    HarmonicOscillator2D s(1);

    integrator integral(p,q,r,s);



    cout<<integral.integrate(20)<<endl;
    return 0;
}

