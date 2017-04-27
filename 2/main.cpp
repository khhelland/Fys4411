#include <iostream>
#include <vmc.h>
#include <time.h>



using namespace std;

int main()
{
    clock_t start,stop;

    start = clock();


    vmc a(2,     // N
          1.7,     // steplength
          1.0,   // alpha
          1,     // omega
          0.3 ); // beta
    a.useImportanceSampling = false;
    a.useInteraction = true;
    a.useJastrow = 1;
    a.run(1e6);

    stop = clock();

    double time = (float)(stop-start)/CLOCKS_PER_SEC;
    cout<<time<<endl;
    return 0;
}

