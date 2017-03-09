#include <iostream>
#include <hartreefock.h>
#include <time.h>
#include <fstream>
//#include <omp.h>

using namespace std;


int main()
{
    clock_t start_total, start_part, end_total, end_part;
    start_total = clock();
    ofstream out("Energies.dat");

    out<<"np,s,w\t"<<"Energy\t"<<"time"<<endl;
    cout<<"np,s,w\t"<<"Energy\t"<<"time"<<endl;

    for(int n = 1; n<5;n++)
    {
        double oldE = 1;
        double newE = 0;
        int j=0;
        while((abs(newE-oldE)>1e-4)&&((n+j)<=12))
        {
            start_part = clock();
            oldE = newE;
            hartreefock hf(n*(n+1),n+j,1);
            hf.run(100,1e-5);
            newE = hf.getenergy();
            end_part = clock();
            double time_part = ((float)(end_part-start_part))/CLOCKS_PER_SEC;
            out<<n*(n+1)<<","<<n+j<<","<<1<<"\t"<<newE<<"\t"<<time_part<<endl;
            cout<<n*(n+1)<<","<<n+j<<","<<1<<"\t"<<newE<<"\t"<<time_part<<endl;
            j++;
        }
    }


    end_total = clock();
    double time = ((float)(end_total-start_total))/CLOCKS_PER_SEC;
    cout<<"total time: "<<time<<endl;
    return 0;
}

