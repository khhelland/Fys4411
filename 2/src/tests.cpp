#include <iostream>
#include <vmc.h>
#include <slatervmc.h>
#include <time.h>
#include <armadillo>


/*
 * Functions for calculating benchmark results and such.
 *
 */

using namespace std;


void manualSteepestDescent(int n,double omega,double alpha, double beta, double gamma,bool jastrow)
{

    slatervmc a( n, //particles
                 0.01,     // steplength
                 omega,   // omega
                 alpha,     // alpha
                 beta,  //beta
                 1);
    a.useJastrow = jastrow;


    a.steepestDescentMan(1e7,gamma);
}

void noIntnoJastrow2(double omega)
{
    vmc a(0.01,
          omega,
          1,
          0.4,
          1);
    a.useJastrow = 0;
    a.useInteraction = 0;
    a.run(1e7,250000);
    a.printResults();
}

void numdiffcomp()
{
    double analyt = 0;
    double numm = 0;

    for(int i = 0; i < 10; i++)
    {
        clock_t start,mid,end;
        start = clock();
        vmc a(0.01,
              1,
              1,
              0.4,
              1);
        a.useJastrow = 0;
        a.useInteraction = 0;
        a.run(1e7,250000);

        mid = clock();

        vmc b(0.01,
              1,
              1,
              0.4,
              1);
        b.useJastrow = 0;
        b.useInteraction = 0;
        b.useNumDiff = 1;
        b.run(1e7,250000);
        end = clock();
        analyt += (double)(mid - start)/CLOCKS_PER_SEC;
        numm += (double)(end - mid)/CLOCKS_PER_SEC;
    }
    analyt/=10;
    numm/=10;
    cout<<"analytical\tnumerical"<<endl<<analyt<<"\t"<<numm<<endl;
}
void noIntnoJastrowN(int n,double omega)
{
    slatervmc a(n,
                0.01,
                omega,
                1,
                0,
                1);
    a.useJastrow = 0;
    a.useInteraction = 0;
    a.run(1e7,250000);
    a.printResults();
}

