#include <harmonicoscillator2d.h>
#include <coulomb_functions.h>
#include <hartreefock.h>









double matrixelement_as(double hw, int ip,int iq, int ir, int is)
{
    HarmonicOscillator2D p(ip);
    HarmonicOscillator2D q(iq);
    HarmonicOscillator2D r(ir);
    HarmonicOscillator2D s(is);

    //spin conservation
    if (!(p.spin + q.spin == r.spin + s.spin))
    {
        return 0;
    }


    double direct = Coulomb_HO(hw,p.n,p.m,q.n,q.m,r.n,r.m,s.n,s.m);
    double exchange = Coulomb_HO(hw,p.n,p.m,q.n,q.m,s.n,s.m,r.n,r.m);
    return direct-exchange;
}


