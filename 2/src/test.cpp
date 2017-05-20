#include <iostream>
#include <slatervmc.h>
#include <time.h>



using namespace std;

int main2()
{


    clock_t start,stop;

    start = clock();


    slatervmc a(2,     // N
          1,     // steplength
          1,     // omega
          0.4 ); // beta
    a.useImportanceSampling = 0;
    a.useInteraction = 1;
    a.useJastrow = 0;

//     cout<<"set up complete"<<endl;


//     a.writeEnergies(1e4,"blocking.dat");
     a.run(1e6, //nCycles
           1);  //Blocksize
     a.printResults();
//     a.steepestDescent(1e5,0.4);

    stop = clock();

    double time = (float)(stop-start)/CLOCKS_PER_SEC;
    cout<<"time [s]: "<<time<<endl;
    return 0;
}

