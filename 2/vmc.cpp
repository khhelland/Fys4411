#include "vmc.h"

#include <armadillo>
#include <random>
#include <iostream>
#include <iomanip>
#include <math.h>


using namespace std;
using namespace arma;

vmc::vmc()
{

}
vmc::vmc(int nParticles, double rSigma, double a, double w):
    nParticles(nParticles)
{
    distribution = normal_distribution<double>(0.0,rSigma);
    uniformDistribution = uniform_real_distribution<double>(0.0,1.0);
    positions = zeros<mat>(nParticles,2);
    positions(0,0) = 1.0;
    alpha = a;
    omega = w;
    generator.seed(5);
    //positions(0,0) = 1;
}

void vmc::run(int nCycles)
{
    double Esum = 0;
    double Esum2 = 0;

    double oldwave = wavefunctionSquared(positions);


    double xstep,ystep,newwave,deltaE,draw,ratio;
    int acceptcount = 0;
    for(int i = 0;i<nCycles;i++)
    {
        /*
        for(int p = 0; p < nParticles; p++)//particles
        {
        */
            int p = 2*(uniformDistribution(generator));
            suggestion = positions;

            xstep = 2*(uniformDistribution(generator) - 0.5);
            ystep = 2*(uniformDistribution(generator) - 0.5);

            suggestion(p,0) += xstep*1.5;
            suggestion(p,1) += ystep*1.5;

            newwave = wavefunctionSquared(suggestion);

            ratio = newwave/oldwave;
            draw = uniformDistribution(generator);

            if(ratio>draw)
            {
                acceptcount++;
                positions = suggestion;
                oldwave = newwave;
            }

            deltaE = localEnergy(positions);

            //save deltaE

            Esum  += deltaE;
            //cout << Esum << endl;
            Esum2 += deltaE*deltaE;
        //}
    }

//    cout << setprecision(15) << Esum << endl;
//    cout << setprecision(15) << Esum2 << endl;

    Esum /=(double)nCycles;
    Esum2/=(double)nCycles;

    cout<< acceptcount/(double)nCycles <<endl;

    //save Esum, Esum2
    cout<<Esum<<","<<Esum2<<","<< Esum2 - Esum*Esum <<endl;

}

double vmc::wavefunctionSquared(mat pos)
{
    return exp(-alpha*omega*accu(pos % pos));
}

double vmc::localEnergy(mat pos)
{
    double xdiff = pos(0,0)-pos(1,0);
    double ydiff = pos(0,1)-pos(1,1);
    return 0.5*(1-alpha*alpha)*omega*omega*accu(pos%pos) + 2.0*alpha*omega + 1.0/sqrt(xdiff*xdiff + ydiff*ydiff);
    //exp(-2*pos(1,0));
    /*

    double l = -0.5*laplacian(pos)/wavefunction(pos) + 0.5*accu(pos%pos);
    //cout << l << endl;
    return l;*/
}

//double vmc::laplacian(mat pos)
//{
//    double sum = 0;
//    double h = 1e-3;
//    double h2 = 1e6;
//    mat change = zeros<mat>(nParticles,2);
//    for(int n = 0;n<nParticles;n++)
//    {
//        for(int d =0;d<2;d++)
//        {
//            change(n,0) = h;
//            sum += wavefunction(pos + change) + wavefunction(pos-change) - 2*wavefunction(pos);
//            change(n,0) = 0;
//        }
//    }
//    return sum*h2;
//}

//double vmc::wavefunction(mat pos)
//{
//    return exp(-accu(pos%pos)/2);
//}
