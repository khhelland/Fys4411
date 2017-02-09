#include <iostream>
#include "harmonicoscillator2d.h"
#include <hermite.h>

using namespace std;

int main2()
{
    HarmonicOscillator2D p(1,1);
    HarmonicOscillator2D q(37,1);

    cout<<p.wavefunction(1,1)<<" "<<q.wavefunction(2,2)<<endl;
    return 0;
}

