#include <iostream>
#include <vmc.h>
#include <time.h>



using namespace std;

int main2()
{


    clock_t start,stop;

    start = clock();


    vmc a(1.5,     // steplength
          0.5,   // alpha
          1,     // omega
          0.4 ); // beta
    a.useImportanceSampling = 0;
    a.useInteraction = 0;
    a.useJastrow = 0;
    a.useNumDiff = 0;



//     a.writeEnergies(1e5,"blocking.dat");
     a.run(1e6, //nCycles
           1);  //Blocksize
     a.printResults();
//     a.steepestDescent(1e5,0.4);

    stop = clock();

    double time = (float)(stop-start)/CLOCKS_PER_SEC;
    cout<<"time [s]: "<<time<<endl;
    return 0;
}

