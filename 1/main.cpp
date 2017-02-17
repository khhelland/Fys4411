#include <iostream>
#include "harmonicoscillator2d.h"
#include <coulomb_functions.h>


using namespace std;

int main()
{
    HarmonicOscillator2D p(109);
    HarmonicOscillator2D q(0);
    HarmonicOscillator2D r(43);
    HarmonicOscillator2D s(16);

//    integrator integral(p,q,r,s);



////    cout<<integral.integrate(20)<<endl;
    double w = 1;
//    int w,e,r,t,y,u,i,o = 1;

    double a = Coulomb_HO(w,p.n,p.m,q.n,q.m,r.n,r.m,s.n,s.m);

    cout<<a<<endl;
    return 0;
}

