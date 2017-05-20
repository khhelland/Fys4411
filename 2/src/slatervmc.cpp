#include "slatervmc.h"

#include <armadillo>
#include <random>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <ho2d.h>


using namespace std;
using namespace arma;

/*
 * do sd and blocking better!!!!!
 * make some funcs inline for opt
 * only non-interacting ground state slater
 *
 * to extend nParticles
 * need new
 * local energy
 * wavefunction and square
 * drift, maybe propratio
 * dirstribute
 * a func -
 * steepest
 *
 */


slatervmc::slatervmc()
{

}
slatervmc::slatervmc(int nParticles, double step, double w, double b):
    nParticles(nParticles),
    stepLength(step),
    omega(w),
    beta(b)
{
    normalDistribution = normal_distribution<double>(0.0,1.0);
    uniformDistribution = uniform_real_distribution<double>(0.0,1.0);
    nOrbitals = nParticles/2;
    distributeParticles();
    setupMatrices();
    olddrift = mat(2,nParticles);
    stepLength2 = stepLength*stepLength;
    generator.seed(1);

}

void slatervmc::distributeParticles()
{
    positions = mat(2,nParticles);
    for(int p = 0;p<nParticles;p++)
    {
        positions(0,p) = p;
        positions(1,p) = p;
    }
}

void slatervmc::setupMatrices()
{
    slaterDown = mat(nOrbitals,nOrbitals);
    slaterUp = mat(nOrbitals,nOrbitals);

    for(int i = 0; i < nOrbitals; i++)
    {
        for(int j = 0; j<nOrbitals; j++)
        {
            slaterDown(i,j) = ho2d(j,omega,positions(0,i),positions(1,i));
            slaterUp(i,j) = ho2d(j,omega,positions(0,i+nOrbitals),positions(1,i+nOrbitals));
        }
    }
    inverseDown = inv(slaterDown);
    inverseUp = inv(slaterUp);
}

void slatervmc::run(int nCycles, int blocksize)
{
    updatePointers();

    updateOld();

    // burn-in
    // for(int i = 0; i<1e5 ; i++) metropolisMove();

    double Esum = 0;
    double Esum2 = 0;
    double Eblock = 0;


    int acceptcount = 0;
    for(int i = 1;i<=nCycles;i++)
    {
        acceptcount += metropolisMove();

        Eblock += (this->*localEnergyPointer)();

        if(i % blocksize == 0)
        {
            Eblock /= (double) blocksize;
            Esum += Eblock;
            Esum2 += Eblock*Eblock;
            Eblock = 0;
        }
    }

    Esum /= (((double)nCycles)/blocksize);
    Esum2 /= (((double)nCycles)/blocksize);

    energy = Esum;
    energySquared = Esum2;
    AcceptanceRatio = acceptcount/(double)nCycles;


}

int slatervmc::metropolisMove()
{
    //randomly choose particle to move along random dimension
    int particle = nParticles*(uniformDistribution(generator));
    int dimension = 2*uniformDistribution(generator);

    suggestion = positions;

    (this->*findSuggestionPointer)(dimension,particle);

    double ratio = (this->*findRatioPointer)(dimension,particle);

    double draw = uniformDistribution(generator);

    if(ratio>draw)
    {
        positions = suggestion;
        olddrift(dimension,particle) = drift;
        updateSlaters(particle);
        return 1;
    }
    else return 0;
}

/*------------------------------------------------------------------------------*/
void slatervmc::findSuggestionUniform(int d,int p)
{
    double step = 2*(uniformDistribution(generator) - 0.5);
    suggestion(d,p) += step*stepLength;
}

void slatervmc::findSuggestionImportanceSamplingwithJastrow(int d, int p)
{
    double step = normalDistribution(generator);
    drift = drifttermwithJastrow(positions,d,p);
    suggestion(d,p) += step*stepLength + drift*stepLength2;
}

void slatervmc::findSuggestionImportanceSamplingnoJastrow(int d, int p)
{
    double step = normalDistribution(generator);
    drift = drifttermnoJastrow(positions,d,p);
    suggestion(d,p) += step*stepLength + drift*stepLength2;
}
/*-------------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------------*/

double slatervmc::findRatio(int d, int p)
{

    double R = 0;
    double xnew = suggestion(0,p);
    double ynew = suggestion(1,p);

    if(p<nOrbitals)
    {
        for(int j = 0; j < nOrbitals; j++)
        {
            R += ho2d(j, omega, xnew, ynew)*inverseDown(j,p);
        }
    }
    else
    {
        for(int j = 0; j <nOrbitals; j++)
        {
            R += ho2d(j, omega, xnew, ynew)*inverseUp(j,p-nOrbitals);
        }
    }
    return R*R;
}

double slatervmc::findRatioJastrow(int d ,int p)
{

    double exponent = 0;
    double rn,ro,av;
    for(int i = 0; i < nParticles; i++)
    {
        if(i != p)
        {
            rn = rDifference(suggestion,i,p);
            ro = rDifference(positions,i,p);
            av = a(i,p);
            exponent += av*(rn/(1+beta*rn) - ro/(1+beta*ro));
        }
    }


    return exp(2*exponent)*findRatio(d,p);
}


double slatervmc::findRatioImportance(int d, int p)
{
    return proposalratio(d,p)*findRatio(d,p);
}

double slatervmc::findRatioImportanceJastrow(int d, int p)
{
    return proposalratio(d,p)*findRatioJastrow(d,p);
}

double slatervmc::wavefunctionSquared(mat pos)
{
    return exp(-omega*accu(pos%pos));
}
double slatervmc::wavefunctionSquaredJastrow(mat pos)
{
    double exponent = 0;
    double r;
    for(int i = 1; i < nParticles; i++)
    {
        for(int j = 0; j < i; j++)
        {
            r = rDifference(pos,i,j);
            exponent += a(i,j)*r/(1+beta*r);
        }
    }
    exponent *= 2;
    return exp(exponent)*wavefunctionSquared(pos);

}

double slatervmc::proposalratio(int d, int p)
{
    return exp(0.5*(drift*drift - olddrift(d,p)*olddrift(d,p))*stepLength2
               + positions(d,p)*(drift + olddrift(d,p))
               - suggestion(d,p)*(drift + olddrift(d,p)));
}

/*-----------------------------------------------------------------------------------*/

double slatervmc::rDifference(mat pos, int p, int q)
{
    double xdiff = pos(0,p) - pos(0,q);
    double ydiff = pos(1,p) - pos(1,q);
    return sqrt(xdiff*xdiff + ydiff*ydiff);
}


/*------------------------------------------------------------------------------------*/
double slatervmc::localEnergy()
{
    return 2.0*omega;
}

double slatervmc::localEnergyJastrow()
{
    double sum = 0;
    double r;
    double b;
    double av;
    for(int i = 1; i < nParticles; i++)
    {
        for(int j = 0; j < i; j++)
        {
            r = rDifference(positions,i,j);
            b = 1.0/(1+beta*r);
            av = a(i,j);
            sum += av*b*b*(-av*b*b - 1.0/r + 2*beta*b + omega*r);

        }
    }
    sum += localEnergy();
    return sum;
}

double slatervmc::localEnergyInteraction()
{
    double sum = 0;
    for(int i = 1; i < nParticles; i++)
    {
        for(int j = 0; j < i; j++)
        {
            sum += 1.0/rDifference(positions,i,j);
        }
    }

    sum += localEnergy();
    return sum;
}

double slatervmc::localEnergyJastrowInteraction()
{
    double sum = 0;
    double r;
    double b;
    double av;
    for(int i = 0; i < nParticles; i++)
    {
        for(int j = i+1; j < nParticles; j++)
        {
            r = rDifference(positions,i,j);
            b = 1.0/(1+beta*r);
            av = a(i,j);
            sum += 1.0/r;
            sum += av*b*b*(-av*b*b - 1.0/r + 2*beta*b + omega*r);

        }
    }

    sum += localEnergy();
    return sum;
}



double slatervmc::wavefunction(mat pos)
{
    return exp(-accu(pos%pos)/2);
}

double slatervmc::wavefunctionJastrow(mat pos)
{
    double r = rDifference(pos,1,0);
    return exp( -accu(pos%pos) + a(1,0)*r/(1+beta*r) );
}

/*---------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
double slatervmc::drifttermnoJastrow(mat pos, int d , int p)
{
    return -omega*pos(d,p);
}

double slatervmc::drifttermwithJastrow(mat pos,int d, int p)
{
    double sum = 0;
    double r;
    for(int q = 0; q < nParticles; q++)
    {
        if(p != q)
        {
            r = rDifference(pos,p,q);
            sum += a(p,q)*(pos(d,p)-pos(d,q))/((1+beta*r)*(1+beta*r)*r);
        }
    }
    return 0.5*(sum -omega*pos(d,p));
}
/*-------------------------------------------------------------------------*/


void slatervmc::writeEnergies(int nCycles, const char* filename)
{
    ofstream out(filename);
    updatePointers();
    updateOld();
    double e;
    for(int i = 0;i<nCycles;i++)
    {
        metropolisMove();
        e = (this->*localEnergyPointer)();
        out<<e<<endl;
    }
    out.close();
}


void slatervmc::steepestDescent(int nCycles, double gamma)
{/*
    if((!useJastrow)||(useNumDiff)||(!useInteraction))
    {
        cout<<"Not implemented for this configuration"<<endl;
    }

    else
    {
        double oldalpha = 0;
        double oldbeta = 0;
        while((abs(oldalpha - alpha) > 1e-5)&&(abs(oldbeta - beta) > 1e-5))
        {
            oldalpha = alpha;
            oldbeta = beta;


            updateOld();


            double Esum = 0;
            double Asum = 0;
            double Bsum = 0;
            double EAsum = 0;
            double EBsum = 0;


            double deltaE,dEda,dEdb;





            for(int i = 0;i<nCycles;i++)
            {
                metropolisMove();
                deltaE = localEnergyJastrowInteraction();
                dEda = alphaDeriv();
                dEdb = betaDeriv();

                Esum += deltaE;
                Asum += dEda;
                Bsum += dEdb;
                EAsum += deltaE*dEda;
                EBsum += deltaE*dEdb;

            }


            Esum /= nCycles;
            Asum /= nCycles;
            Bsum /= nCycles;
            EAsum /= nCycles;
            EBsum /= nCycles;


            energy = Esum;
            alpha -= gamma*2*(EAsum -Esum*Asum);
            beta -= gamma*2*(EBsum -Esum*Bsum);
            cout<<alpha<<", "<<beta<<endl;


        }
    }*/
}

double slatervmc::betaDeriv()
{
    double sum = 0;
    double r;
    double s;
    for(int i = 1; i<nParticles; i++)
    {
        for(int j = 0; j<i; j++)
        {
            r = rDifference(positions,i,j);
            s = r/(1+beta*r);
            sum -= a(i,j)*s*s;
        }
    }

    return sum;

}

void slatervmc::printResults()
{
    cout<< "Acceptanceratio: "<<AcceptanceRatio <<endl;

    //save Esum, Esum2
    cout<<"E, E^2, sigma" <<endl;
    cout<<energy<<","<<energySquared<<","<< sqrt(abs(energySquared - energy*energy)) <<endl;
}

void slatervmc::updateSlaters(int p)
{
    double xnew = positions(0,p);
    double ynew = positions(1,p);

    if(p<nOrbitals)
    {

        for(int j = 0; j<nOrbitals; j++)
        {
            slaterDown(p,j) = ho2d(j,omega, xnew, ynew);
        }
        //Consider inserting efficient algorithm here
        inverseDown = inv(slaterDown);
    }
    else
    {
        for(int j = 0; j<nOrbitals; j++)
        {
            slaterUp(p-nOrbitals,j) = ho2d(j,omega, xnew, ynew);
        }
        //Consider inserting efficient algorithm here
        inverseUp = inv(slaterUp);

    }
}

void slatervmc::updateOld()
{
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

    if(useJastrow) oldwavesquared = wavefunctionSquaredJastrow(positions);
    else oldwavesquared = wavefunctionSquared(positions);
}

void slatervmc::updatePointers()
{
    // Make sure function pointers match current booleans
    // 8 cases
    /*------------------------------------------------------------------------------*/


    if ((useImportanceSampling)&&(useInteraction)&&(useJastrow))
    {
        findSuggestionPointer = &slatervmc::findSuggestionImportanceSamplingwithJastrow;
        findRatioPointer = &slatervmc::findRatioImportanceJastrow;
        localEnergyPointer = &slatervmc::localEnergyJastrowInteraction;
    }

    else if ((useImportanceSampling)&&(useInteraction)&&(!useJastrow) )
    {
        findSuggestionPointer = &slatervmc::findSuggestionImportanceSamplingnoJastrow;
        findRatioPointer = &slatervmc::findRatioImportance;
        localEnergyPointer = &slatervmc::localEnergyInteraction;
    }

    else if ((useImportanceSampling)&&(!useInteraction)&&(useJastrow) )
    {
        findSuggestionPointer = &slatervmc::findSuggestionImportanceSamplingwithJastrow;
        findRatioPointer = &slatervmc::findRatioImportanceJastrow;
        localEnergyPointer = &slatervmc::localEnergyJastrow;
    }

    else if ((!useImportanceSampling)&&(useInteraction)&&(useJastrow) )
    {
        findSuggestionPointer = &slatervmc::findSuggestionUniform;
        findRatioPointer = &slatervmc::findRatioJastrow;
        localEnergyPointer = &slatervmc::localEnergyJastrowInteraction;
    }


    else if ((!useImportanceSampling)&&(!useInteraction)&&(useJastrow) )
    {
        findSuggestionPointer = &slatervmc::findSuggestionUniform;
        findRatioPointer = &slatervmc::findRatioJastrow;
        localEnergyPointer = &slatervmc::localEnergyJastrow;
    }

    else if ((!useImportanceSampling)&&(useInteraction)&&(!useJastrow) )
    {
        findSuggestionPointer = &slatervmc::findSuggestionUniform;
        findRatioPointer = &slatervmc::findRatio;
        localEnergyPointer = &slatervmc::localEnergyInteraction;
    }

    else if ((useImportanceSampling)&&(!useInteraction)&&(!useJastrow) )
    {
        findSuggestionPointer = &slatervmc::findSuggestionImportanceSamplingnoJastrow;
        findRatioPointer = &slatervmc::findRatioImportance;
        localEnergyPointer = &slatervmc::localEnergy;
    }

    else if ((!useImportanceSampling)&&(!useInteraction)&&(!useJastrow) )
    {
        findSuggestionPointer = &slatervmc::findSuggestionUniform;
        findRatioPointer = &slatervmc::findRatio;
        localEnergyPointer = &slatervmc::localEnergy;
    }

    /*-------------------------------------------------------------------------------*/

}
