#include <iostream>
#include <hartreefock.h>


using namespace std;


int main()
{
    hartreefock sys(6,3,1);


    sys.run(100,1e-5);
//    sys.fockmatrix.print();
    sys.print_sp_energy();
    cout<<endl<<sys.getenergy()<<endl;

    return 0;
}

