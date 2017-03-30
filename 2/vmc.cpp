#include "vmc.h"

#include <armadillo>
#include <random>
#include <iostream>
#include <math.h>

using namespace std;
using namespace arma;

vmc::vmc()
{

}
vmc::vmc(int nParticles, double rSigma):
    nParticles(nParticles)
{
    distribution = normal_distribution<double>(0.0,rSigma);
    uniformDistribution = uniform_real_distribution<double>(0.0,1.0);
    positions = zeros<mat>(nParticles,2);
}

void vmc::run(int nCycles)
{
    double Esum = 0;
    double Esum2 = 0;

    double oldwave = wavefunctionSquared(positions);


    double xstep,ystep,newwave,deltaE,draw,ratio;
    for(int i = 1;i<nCycles;i++)
    {
        for(int p = 0; p < nParticles; p++)//particles
        {
            suggestion = positions;

            xstep = distribution(generator);
            ystep = distribution(generator);

            suggestion(p,0) += xstep;
            suggestion(p,1) += ystep;

            newwave = wavefunctionSquared(suggestion);

            ratio = newwave/oldwave;
            draw = uniformDistribution(generator);
            if(ratio>draw)
            {
                positions = suggestion;
                oldwave = newwave;
            }

            deltaE = localEnergy(positions);

            //save deltaE

            Esum += deltaE;
            Esum2 += deltaE*deltaE;
        }
    }
    Esum/=nParticles*nCycles;
    Esum2/=nParticles*nCycles;
    //save Esum, Esum2
    cout <<Esum<<","<<Esum2<<","<<Esum2 - Esum*Esum<<endl;
}

double vmc::wavefunctionSquared(mat pos)
{
    return exp(-accu(pos % pos));
}

double vmc::localEnergy(arma::mat pos)
{
    return exp(-2*pos(0,0));
}
