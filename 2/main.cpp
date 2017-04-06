#include <iostream>
#include <vmc.h>
#include <time.h>



using namespace std;

int main()
{
    clock_t start,stop;

    start = clock();


    //    N   sigma  alpha   omega
    vmc a(2,  2,     0.6867,   1);
    a.run(1e8);

    stop = clock();

    double time = (float)(stop-start)/CLOCKS_PER_SEC;
    cout<<time<<endl;
    return 0;
}

