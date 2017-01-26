#include <iostream>
#include "harmonicoscillator2d.h"
#include <hermite.h>

using namespace std;

int main()
{
    HarmonicOscillator2D HO;
    int i = HO.getindex(1,2,-1);

    Hermite poly0(0);
    Hermite poly2(2);


    cout <<poly0.evaluate(2)<<" "<< poly2.evaluate(2)<< endl;
    return 0;
}

