#include <iostream>
#include <slatervmc.h>
#include <time.h>



using namespace std;

int main2()
{


    clock_t start,stop;

    start = clock();


    slatervmc a(12,     // N
          0.01,         // steplength
          1,            // omega
          2,            //alpha
          0.1 );        // beta
    a.useImportanceSampling = 1;
    a.useInteraction = 1;
    a.useJastrow = 1;

//     cout<<"set up complete"<<endl;


//     a.writeEnergies(1e4,"blocking.dat");
//     a.run(1e6, //nCycles
//           1);  //Blocksize
//     a.printResults();
     a.steepestDescent(1e5,0.4);

    stop = clock();

    double time = (float)(stop-start)/CLOCKS_PER_SEC;
    cout<<"time [s]: "<<time<<endl;
    return 0;
}

