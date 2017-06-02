#include <iostream>
#include <vmc.h>
#include <time.h>
#include <tests.h>


using namespace std;

int main()
{


    clock_t start,stop;

    start = clock();
    numdiffcomp();
    //noIntnoJastrow2(1);
//    manualSteepestDescent(6,0.01,1,0.5,0.1,1);
//    vmc a(0.01,     // steplength
//          1,   // omega
//          1,     // alpha
//          0.35,  //beta
//          2);
//    a.useImportanceSampling = 1;
//    a.useInteraction = 1;
//    a.useJastrow = 1;
//    a.useNumDiff = 0;

//    cout<<"set up complete"<<endl;

// //    a.writeEnergies(1e5,"blocking.dat");
//     a.run(1e7, //nCycles
//           250000);  //Blocksize
//     a.printResults();
//     cout<<a.mean_distance<<endl;
//    a.steepestDescent(1e5,1);
//    a.run(1e6,5000);
//    a.printResults();

    stop = clock();

    double time = (float)(stop-start)/CLOCKS_PER_SEC;
    cout<<"time [s]: "<<time<<endl;
    return 0;
}

