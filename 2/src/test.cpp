#include <iostream>
#include <slatervmc.h>
#include <time.h>
#include "ho2d.h"
#include "tests.h"

using namespace std;

int main()
{


    clock_t start,stop;

    start = clock();


    slatervmc a(2,     // N
          0.01,         // steplength
          1,            // omega
          1,            //alpha
          0.4,        // beta
          3);         //seed
    a.useImportanceSampling = 1;
    a.useInteraction = 1;
    a.useJastrow = 1;

//     cout<<"set up complete"<<endl;


//     a.writeEnergies(1e7,"e2pIsJ.dat");
     a.run(1e7, //nCycles
           2e5);  //Blocksize
     a.printResults();
//     a.steepestDescentAuto(1e5,1);

//     cout<<ho2ddw(7,1,1,1)<<endl;
//    alphaderivcomp();
    stop = clock();

    double time = (float)(stop-start)/CLOCKS_PER_SEC;
    cout<<"time [s]: "<<time<<endl;
    return 0;
}

