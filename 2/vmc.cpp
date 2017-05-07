#include "vmc.h"

#include <armadillo>
#include <random>
#include <iostream>
#include <iomanip>
#include <math.h>


using namespace std;
using namespace arma;

/*
 * Tenk paa aa legge inn burn in
 * member function pointers
 * eller boost::function
 *
 */


vmc::vmc()
{

}
vmc::vmc(int nParticles, double step, double alph, double w, double b):
    nParticles(nParticles),
    stepLength(step),
    alpha(alph),
    omega(w),
    beta(b)
{
    normalDistribution = normal_distribution<double>(0.0,1.0);
    uniformDistribution = uniform_real_distribution<double>(0.0,1.0);
    distributeParticles();
    stepLength2 = stepLength*stepLength;
    generator.seed(20);

}

void vmc::distributeParticles()
{
    positions = mat(2,nParticles);
    for(int p = 0;p<nParticles;p++)
    {
        positions(0,p) = p;
        positions(1,p) = p;
    }
}

void vmc::run(int nCycles, int blocksize)
{
    double Esum = 0;
    double Esum2 = 0;
    double Eblock = 0;

    double oldwave = wavefunctionSquared(positions);
    mat olddrift(2,nParticles);
    if(useImportanceSampling)
    {
        for(int i = 0;i<nParticles;i++)
        {
            for(int j=0;j<2;j++)
            {
                olddrift(j,i) = driftterm(positions,i,j);
            }
        }
    }


    double step,newwave,deltaE,draw,ratio, proposalratio;
    double drift = 0;
    int particle,dimension;



    int acceptcount = 0;
    for(int i = 0;i<nCycles;i++)
    {
        //randomly choose particle to move along random dimension
        particle = nParticles*(uniformDistribution(generator));
        dimension = 2*uniformDistribution(generator);

        suggestion = positions;

        if(useImportanceSampling)
        {
            //cout<<"Importance Sampling"<<endl;
            step = normalDistribution(generator);
        }
        else
        {
            step = 2*(uniformDistribution(generator)-0.5);
        }

        suggestion(dimension,particle) += step*stepLength;


        ratio = 1;
        if (useImportanceSampling)
        {
            drift = driftterm(positions,particle,dimension);
            suggestion(dimension,particle) += stepLength2*drift;
            proposalratio = proposalDensity(positions(dimension,particle),  suggestion(dimension,particle), olddrift(dimension,particle)) /
                            proposalDensity(suggestion(dimension,particle), positions(dimension,particle),  drift) ;
            ratio = proposalratio;
        }

        //Metropolis algorithm:
        newwave = wavefunctionSquared(suggestion);

        ratio  *= newwave/oldwave;

        draw = uniformDistribution(generator);

        if(ratio>draw)
        {
            acceptcount++;
            positions = suggestion;
            oldwave = newwave;       
            olddrift(dimension,particle) = drift;
        }

        deltaE = localEnergy(positions);

        //cout<<deltaE<<endl;

        //save deltaE


        Eblock += deltaE;

        if(i % blocksize == 0)
        {
            Esum += Eblock;
            Esum2 += Eblock*Eblock;
            Eblock = 0;
        }
    }

    //    cout << setprecision(15) << Esum << endl;
    //    cout << setprecision(15) << Esum2 << endl;

    Esum /=nCycles;
    Esum2/=nCycles;

    cout<< "Acceptanceratio: "<<acceptcount/(double)nCycles <<endl;

    //save Esum, Esum2
    cout<<"E, E^2, sigma" <<endl;
    cout<<Esum<<","<<Esum2<<","<< sqrt(Esum2 - Esum*Esum) <<endl;

}

double vmc::rDifference(mat pos, int p, int q)
{
    double xdiff = pos(0,p) - pos(0,q);
    double ydiff = pos(1,p) - pos(1,q);
    return sqrt(xdiff*xdiff + ydiff*ydiff);
}

double vmc::wavefunctionSquared(mat pos)
{
    double w = 1;

    if (useJastrow)
    {
        double r;
        for(int i=1;i<nParticles;i++)
        {
            for(int j=0; j<i;j++)
            {
                r = rDifference(pos,i,j);
                w *=exp(2*a*r/(1+beta*r));
            }
        }
    }

    return w*exp(-alpha*omega*accu(pos % pos));
}

double vmc::localEnergy(mat pos)
{
    double sum = 0;
    sum  += 0.5*(1-alpha*alpha)*omega*omega*accu(pos%pos) + 2.0*alpha*omega;
    //cout << sum << endl;
//    if(useInteraction)
//    {
//        sum += 1/rDifference(pos,0,1);
//    }
//    if(useJastrow)
//    {
//        double r = rDifference(pos,0,1);
//        double b = 1/(1+beta*r);
//        sum += -a*b*b*(a*b*b + 1.0/r - 2*beta*b - alpha*omega*r) ;
//    }
    if(useInteraction&&useJastrow)
    {
        //cout <<"both"<<", "<<beta<<endl;
        double r;
        double b;
        for(int i = 0; i < nParticles; i++)
        {
            for(int j = i+1; j < nParticles; j++)
            {
                r = rDifference(pos,i,j);
                b = 1.0/(1+beta*r);
                sum += 1.0/r;
                sum += a*b*b*(-a*b*b - 1.0/r + 2*beta*b + alpha*omega*r);

            }
        }
    }

    else if(useInteraction&&(!useJastrow))
    {
        //cout<<"int"<<endl;
        for(int i = 1; i < nParticles; i++)
        {
            for(int j = 0; j < i; j++)
            {
                //cout<<"particles"<<i<<","<<j<<", " << 1.0/rDifference(i,j)<<endl;
                sum += 1.0/rDifference(pos,i,j);


            }
        }
    }

    else if((!useInteraction)&&useJastrow)
    {
        //cout<<"jastrow"<<endl;
        double r;
        double b;
        for(int i = 1; i < nParticles; i++)
        {
            for(int j = 0; j < i; j++)
            {
                r = rDifference(pos,i,j);
                b = 1.0/(1+beta*r);
                sum += a*b*b*(-a*b*b - 1.0/r + 2*beta*b + alpha*omega*r);

            }
        }
    }
    return sum;
}

double vmc::driftterm(mat pos, int p, int dim)  //(nabla psi)/2*psi
{
    double d = -2*alpha*omega*pos(dim,p); //slater
    if(useJastrow)
    {
        double r;
        for(int q = 0; q < nParticles; q++)
        {
            if(p != q)
            {
                r = rDifference(pos,p,q);
                d += a*(pos(dim,p)-pos(dim,q))/((1+beta*r)*(1+beta*r)*r);
            }
        }
    }

    return 0.5*d;
}

double vmc::proposalDensity(double rn, double ro, double drift)
{
    return exp(-(rn-ro-stepLength2*drift)*(rn-ro-stepLength2*drift)/(2*stepLength2));
}


double vmc::laplacian(mat pos)
{
    double sum = 0;
    double h = 1e-3;
    double h2 = 1e6;
    mat change = zeros<mat>(nParticles,2);
    for(int n = 0;n<nParticles;n++)
    {
        for(int d =0;d<2;d++)
        {
            change(n,d) = h;
            sum += wavefunction(pos + change) + wavefunction(pos-change);
            change(n,d) = 0;
        }
    }
    sum -= 4*nParticles*wavefunction(pos);
    return sum*h2;
}

double vmc::wavefunction(mat pos)
{
    return exp(-accu(pos%pos)/2);
}
