#include <iostream>
#include "harmonicoscillator2d.h"
#include <hartreefock.h>


using namespace std;


int main()
{
    hartreefock sys(2,5,1);

//    sys.ref_energies.print();

//    cout<<endl;
//    cout<<endl;
//    sys.densitymatrix.print();

//    cout<<endl;
//    cout<<endl;
//    sys.fockmatrix.print();

//    sys.updatefockmatrix();

//    cout<<endl;
//    cout<<endl;
//    sys.fockmatrix.print();

    sys.run(100,1e-5);
    cout<<sys.getenergy()<<endl;

    return 0;
}

