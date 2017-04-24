#include <iostream>
#include <vmc.h>
#include <time.h>



using namespace std;

int main()
{
    clock_t start,stop;

    start = clock();


    //    N   steplength  alpha   omega
    vmc a(2,  0.05,     0.7,   1);
    a.useImportanceSampling = false;
    a.useInteraction = false;
    a.useJastrow = false;
    a.run(1e6);

    stop = clock();

    double time = (float)(stop-start)/CLOCKS_PER_SEC;
    cout<<time<<endl;
    return 0;
}

