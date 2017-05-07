#include "mc.h"
#include <iostream>


using namespace std;

mc::mc()
{
}

int mc::neg(int n)
{
    return -n;
}

int mc::pos(int n)
{
    return n;
}

void mc::val(int n)
{
    int (mc::*p)(int);
    if(n%2==0) p = &mc::pos;
    else p = &mc::neg;
    cout<<(this->*p)(n)<<endl;
}
