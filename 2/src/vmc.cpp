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
 * Think about saving oldwave again.
 *Egen run for steepest descent?
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
    olddrift = mat(2,nParticles);
    stepLength2 = stepLength*stepLength;
    generator.seed(1);

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
    // First set up function pointers
    // 16 cases
    /*------------------------------------------------------------------------------*/
    void (vmc::*findSuggestionPointer)(int,int);
    double (vmc::*findRatioPointer)(int,int);
    double (vmc::*localEnergyPointer)();

    if ((useImportanceSampling)&&(useInteraction)&&(useJastrow)&&(useNumDiff))
    {
        findSuggestionPointer = &vmc::findSuggestionImportanceSamplingwithJastrow;
        findRatioPointer = &vmc::findRatioImportanceJastrow;
        localEnergyPointer = &vmc::localEnergyJastrowInteractionNumdiff;
    }

    else if ((useImportanceSampling)&&(useInteraction)&&(useJastrow)&&(!useNumDiff))
    {
        findSuggestionPointer = &vmc::findSuggestionImportanceSamplingwithJastrow;
        findRatioPointer = &vmc::findRatioImportanceJastrow;
        localEnergyPointer = &vmc::localEnergyJastrowInteraction;
    }

    else if ((useImportanceSampling)&&(useInteraction)&&(!useJastrow)&&(useNumDiff))
    {
        findSuggestionPointer = &vmc::findSuggestionImportanceSamplingnoJastrow;
        findRatioPointer = &vmc::findRatioImportance;
        localEnergyPointer = &vmc::localEnergyInteractionNumdiff;
    }

    else if ((useImportanceSampling)&&(!useInteraction)&&(useJastrow)&&(useNumDiff))
    {
        findSuggestionPointer = &vmc::findSuggestionImportanceSamplingwithJastrow;
        findRatioPointer = &vmc::findRatioImportanceJastrow;
        localEnergyPointer = &vmc::localEnergyJastrowNumdiff;
    }

    else if ((!useImportanceSampling)&&(useInteraction)&&(useJastrow)&&(useNumDiff))
    {
        findSuggestionPointer = &vmc::findSuggestionUniform;
        findRatioPointer = &vmc::findRatioJastrow;
        localEnergyPointer = &vmc::localEnergyJastrowInteractionNumdiff;
    }

    else if ((useImportanceSampling)&&(useInteraction)&&(!useJastrow)&&(!useNumDiff))
    {
        findSuggestionPointer = &vmc::findSuggestionImportanceSamplingnoJastrow;
        findRatioPointer = &vmc::findRatioImportance;
        localEnergyPointer = &vmc::localEnergyInteraction;
    }

    else if ((useImportanceSampling)&&(!useInteraction)&&(useJastrow)&&(!useNumDiff))
    {
        findSuggestionPointer = &vmc::findSuggestionImportanceSamplingwithJastrow;
        findRatioPointer = &vmc::findRatioImportanceJastrow;
        localEnergyPointer = &vmc::localEnergyJastrow;
    }

    else if ((!useImportanceSampling)&&(useInteraction)&&(useJastrow)&&(!useNumDiff))
    {
        findSuggestionPointer = &vmc::findSuggestionUniform;
        findRatioPointer = &vmc::findRatioJastrow;
        localEnergyPointer = &vmc::localEnergyJastrowInteraction;
    }

    else if ((useImportanceSampling)&&(!useInteraction)&&(!useJastrow)&&(useNumDiff))
    {
        findSuggestionPointer = &vmc::findSuggestionImportanceSamplingnoJastrow;
        findRatioPointer = &vmc::findRatioImportance;
        localEnergyPointer = &vmc::localEnergyNumdiff;
    }


    else if ((!useImportanceSampling)&&(useInteraction)&&(!useJastrow)&&(useNumDiff))
    {
        findSuggestionPointer = &vmc::findSuggestionUniform;
        findRatioPointer = &vmc::findRatio;
        localEnergyPointer = &vmc::localEnergyInteractionNumdiff;
    }

    else if ((!useImportanceSampling)&&(!useInteraction)&&(useJastrow)&&(useNumDiff))
    {
        findSuggestionPointer = &vmc::findSuggestionUniform;
        findRatioPointer = &vmc::findRatioJastrow;
        localEnergyPointer = &vmc::localEnergyJastrowNumdiff;
    }

    else if ((!useImportanceSampling)&&(!useInteraction)&&(!useJastrow)&&(useNumDiff))
    {
        findSuggestionPointer = &vmc::findSuggestionUniform;
        findRatioPointer = &vmc::findRatio;
        localEnergyPointer = &vmc::localEnergyNumdiff;
    }

    else if ((!useImportanceSampling)&&(!useInteraction)&&(useJastrow)&&(!useNumDiff))
    {
        findSuggestionPointer = &vmc::findSuggestionUniform;
        findRatioPointer = &vmc::findRatioJastrow;
        localEnergyPointer = &vmc::localEnergyJastrow;
    }

    else if ((!useImportanceSampling)&&(useInteraction)&&(!useJastrow)&&(!useNumDiff))
    {
        findSuggestionPointer = &vmc::findSuggestionUniform;
        findRatioPointer = &vmc::findRatio;
        localEnergyPointer = &vmc::localEnergyInteraction;
    }

    else if ((useImportanceSampling)&&(!useInteraction)&&(!useJastrow)&&(!useNumDiff))
    {
        findSuggestionPointer = &vmc::findSuggestionImportanceSamplingnoJastrow;
        findRatioPointer = &vmc::findRatioImportance;
        localEnergyPointer = &vmc::localEnergy;
    }

    else if ((!useImportanceSampling)&&(!useInteraction)&&(!useJastrow)&&(!useNumDiff))
    {
        findSuggestionPointer = &vmc::findSuggestionUniform;
        findRatioPointer = &vmc::findRatio;
        localEnergyPointer = &vmc::localEnergy;
    }

    /*-------------------------------------------------------------------------------*/


    double Esum = 0;
    double Esum2 = 0;
    double Eblock = 0;


    mat olddrift(2,nParticles);
    if(useImportanceSampling)
    {
        for(int p = 0; p<nParticles; p++)
        {
            for(int d=0;d<2;d++)
            {
                if(useJastrow)
                    olddrift(d,p) = drifttermwithJastrow(positions,d,p);
                else
                    olddrift(d,p) = drifttermnoJastrow(positions,d,p);
            }
        }
    }


    double deltaE,draw,ratio;
    int particle,dimension;



    int acceptcount = 0;
    for(int i = 0;i<nCycles;i++)
    {
        //randomly choose particle to move along random dimension
        particle = nParticles*(uniformDistribution(generator));
        dimension = 2*uniformDistribution(generator);

        suggestion = positions;
        // cout<<"entered loop"<<endl;

        (this->*findSuggestionPointer)(dimension,particle);

        // cout<<"found suggestion"<<endl;

        ratio = (this->*findRatioPointer)(dimension,particle);

        // cout<<"found ratio"<<endl;


        draw = uniformDistribution(generator);

        if(ratio>draw)
        {
            // (this->*updateSystem)();
            acceptcount++;
            positions = suggestion;                 
            olddrift(dimension,particle) = drift;
        }

        deltaE = (this->*localEnergyPointer)();

        //cout<<deltaE<<endl;


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

    energy = Esum;
    energySquared = Esum2;
    AcceptanceRatio = acceptcount/(double)nCycles;


}

/*------------------------------------------------------------------------------*/
void vmc::findSuggestionUniform(int d,int p)
{
    double step = 2*(uniformDistribution(generator) - 0.5);
    suggestion(d,p) += step*stepLength;
}

void vmc::findSuggestionImportanceSamplingwithJastrow(int d, int p)
{
    double step = normalDistribution(generator);
    drift = drifttermwithJastrow(positions,d,p);
    suggestion(d,p) += step*stepLength + drift*stepLength2;
}

void vmc::findSuggestionImportanceSamplingnoJastrow(int d, int p)
{
    double step = normalDistribution(generator);
    drift = drifttermnoJastrow(positions,d,p);
    suggestion(d,p) += step*stepLength + drift*stepLength2;
}
/*-------------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------------*/
double vmc::findRatioJastrow(int d ,int p)
{
    double exponent = 0;
    double rnew;
    double rold;
    for(int i = 1; i < nParticles; i++)
    {
        for(int j = 0; j < i; j++)
        {
            rnew = rDifference(suggestion,i,j);
            rold = rDifference(positions,i,j);
            exponent += (rnew/(1+beta*rnew) - rold/(1+beta*rold));
        }
    }
    exponent *= 2*a;
    exponent += alpha*omega*(-accu(suggestion%suggestion) + accu(positions%positions));
    return exp(exponent);
}

double vmc::findRatio(int d, int p)
{
    return exp(alpha*omega*(-accu(suggestion%suggestion) + accu(positions%positions)));
}

double vmc::findRatioImportance(int d, int p)
{
    double propexponent = 0.5*(drift*drift - olddrift(d,p)*olddrift(d,p))*stepLength2
                          + positions(d,p)*(drift + olddrift(d,p))
                          - suggestion(d,p)*(drift + olddrift(d,p));
    double waveexponent = alpha*omega*(-accu(suggestion%suggestion) + accu(positions%positions));
    return exp(propexponent + waveexponent);
}

double vmc::findRatioImportanceJastrow(int d, int p)
{
    double propexponent = 0.5*(drift*drift - olddrift(d,p)*olddrift(d,p))*stepLength2
                          + positions(d,p)*(drift + olddrift(d,p))
                          - suggestion(d,p)*(drift + olddrift(d,p));
    double waveexponent = 0;
    double rnew;
    double rold;
    for(int k = 1; k < nParticles; k++)
    {
        for(int q = 0; q < k; q++)
        {
            rnew = rDifference(suggestion,k,q);
            rold = rDifference(positions,k,q);
            waveexponent += (rnew/(1+beta*rnew) - rold/(1+beta*rold));
        }
    }
    waveexponent *= 2*a;
    waveexponent += alpha*omega*(-accu(suggestion%suggestion) + accu(positions%positions));

    return exp(propexponent + waveexponent);

}

/*-----------------------------------------------------------------------------------*/

double vmc::rDifference(mat pos, int p, int q)
{
    double xdiff = pos(0,p) - pos(0,q);
    double ydiff = pos(1,p) - pos(1,q);
    return sqrt(xdiff*xdiff + ydiff*ydiff);
}


/*----------------------------------------------------------------*/
double vmc::localEnergy()
{
    return 0.5*(1-alpha*alpha)*omega*omega*accu(positions%positions) + 2.0*alpha*omega;
}

double vmc::localEnergyJastrow()
{
    double sum = 0;
    double r;
    double b;
    for(int i = 1; i < nParticles; i++)
    {
        for(int j = 0; j < i; j++)
        {
            r = rDifference(positions,i,j);
            b = 1.0/(1+beta*r);
            sum += b*b*(-a*b*b - 1.0/r + 2*beta*b + alpha*omega*r);

        }
    }
    sum *= a;
    sum += 0.5*(1-alpha*alpha)*omega*omega*accu(positions%positions) + 2.0*alpha*omega;
    return sum;
}

double vmc::localEnergyInteraction()
{
    double sum = 0;
    for(int i = 1; i < nParticles; i++)
    {
        for(int j = 0; j < i; j++)
        {
            //cout<<"particles"<<i<<","<<j<<", " << 1.0/rDifference(i,j)<<endl;
            sum += 1.0/rDifference(positions,i,j);
        }
    }

    sum += 0.5*(1-alpha*alpha)*omega*omega*accu(positions%positions) + 2.0*alpha*omega;
    return sum;
}

double vmc::localEnergyJastrowInteraction()
{
    double sum = 0;
    double r;
    double b;
    for(int i = 0; i < nParticles; i++)
    {
        for(int j = i+1; j < nParticles; j++)
        {
            r = rDifference(positions,i,j);
            b = 1.0/(1+beta*r);
            sum += 1.0/r;
            sum += a*b*b*(-a*b*b - 1.0/r + 2*beta*b + alpha*omega*r);

        }
    }

    sum += 0.5*(1-alpha*alpha)*omega*omega*accu(positions%positions) + 2.0*alpha*omega;
    return sum;
}

double vmc::localEnergyNumdiff()
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
            sum += wavefunction(positions + change) + wavefunction(positions-change);
            change(n,d) = 0;
        }
    }
    sum /= wavefunction(positions);
    sum -= 4*nParticles;
    return 0.5*(omega*omega*accu(positions%positions)-h2*sum);
}

double vmc::localEnergyJastrowNumdiff()
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
            sum += wavefunctionJastrow(positions + change) + wavefunctionJastrow(positions-change);
            change(n,d) = 0;
        }
    }
    sum /= wavefunctionJastrow(positions);
    sum -= 4*nParticles;
    return 0.5*(omega*omega*accu(positions%positions)-h2*sum);

}

double vmc::localEnergyInteractionNumdiff()
{
    double k = localEnergyNumdiff();
    return k + 1/rDifference(positions,0,1);

}

double vmc::localEnergyJastrowInteractionNumdiff()
{
    double k = localEnergyJastrowNumdiff();
    return k + 1/rDifference(positions,1,0);
}


double vmc::wavefunction(mat pos)
{
    return exp(-accu(pos%pos)/2);
}

double vmc::wavefunctionJastrow(mat pos)
{
    double r = rDifference(pos,1,0);
    return exp( -accu(pos%pos) + a*r/(1+beta*r) );
}

/*---------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
double vmc::drifttermnoJastrow(mat pos, int d , int p)
{
    return -alpha*omega*pos(d,p);
}

double vmc::drifttermwithJastrow(mat pos,int d, int p)
{
    double sum = 0;
    double r;
    for(int q = 0; q < nParticles; q++)
    {
        if(p != q)
        {
            r = rDifference(pos,p,q);
            sum += (pos(d,p)-pos(d,q))/((1+beta*r)*(1+beta*r)*r);
        }
    }
    return 0.5*(a*sum - alpha*omega*pos(d,p));
}
/*-------------------------------------------------------------------------*/


int vmc::blocksize(int nCycles)
{

    // gives 341 for step 10
    // but 3401 for step 100
    // weirdS
    int block = 1;
    int step = 100;
    run(nCycles,block);
    double var =  energySquared - energy*energy;
    double oldvar = 0;
    //course grain
    while(var>oldvar)
    {
        oldvar = var;
        block += step;
        run(nCycles,block);
        var = energySquared - energy*energy;
    }
    //needs finer graining.


    return block;


}


void vmc::printResults()
{
    cout<< "Acceptanceratio: "<<AcceptanceRatio <<endl;

    //save Esum, Esum2
    cout<<"E, E^2, sigma" <<endl;
    cout<<energy<<","<<energySquared<<","<< sqrt(abs(energySquared - energy*energy)) <<endl;
}