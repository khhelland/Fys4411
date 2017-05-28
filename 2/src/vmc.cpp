#include "vmc.h"

#include <armadillo>
#include <random>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>


using namespace std;
using namespace arma;

/*
 * do sd and blocking better!!!!!
 * make some funcs inline for opt
 *
 */


vmc::vmc()
{

}
vmc::vmc(double step, double w, double a, double b, int seed):
    stepLength(step),
    omega(w),
    alpha(a),
    beta(b)
{
    normalDistribution = normal_distribution<double>(0.0,1.0);
    uniformDistribution = uniform_real_distribution<double>(0.0,1.0);
    distributeParticles();
    olddrift = mat(2,nParticles);
    stepLength2 = stepLength*stepLength;
    generator.seed(seed);

}

void vmc::distributeParticles()
{
    positions = mat(2,nParticles,fill::randu);
//    for(int p = 0;p<nParticles;p++)
//    {
//        positions(0,p) = p;
//        positions(1,p) = p;
//    }
}

void vmc::run(int nCycles, int blocksize)
{
    updatePointers();
    updateOld();

    // burn-in
    for(int i = 0; i<1e5 ; i++) metropolisMove();

    double Esum = 0;
    double Esum2 = 0;
    double Eblock = 0;
    double distance_sum = 0;

    int acceptcount = 0;
    for(int i = 1;i<=nCycles;i++)
    {
        acceptcount += metropolisMove();

        Eblock += (this->*localEnergyPointer)();

        distance_sum += rDifference(positions,1,2);

        if(i % blocksize == 0)
        {
            Eblock /= (double) blocksize;
            Esum += Eblock;
            Esum2 += Eblock*Eblock;
            Eblock = 0;
        }
    }

    double  nBlocks = (double)nCycles/blocksize;
    Esum /= nBlocks;
    Esum2 /= nBlocks;

    distance_sum /= nCycles;

    energy = Esum;
    energySquared = Esum2;
    meanVar = (Esum2 - Esum*Esum)/(nBlocks-1);

    mean_distance = distance_sum;


    AcceptanceRatio = acceptcount/(double)nCycles;


}

int vmc::metropolisMove()
{
    //randomly choose particle to move along random dimension
    int particle = nParticles*(uniformDistribution(generator));
    int dimension = 2*uniformDistribution(generator);

    suggestion = positions;
    // cout<<"entered loop"<<endl;

    (this->*findSuggestionPointer)(dimension,particle);

    // cout<<"found suggestion"<<endl;

    double ratio = (this->*findRatioPointer)(dimension,particle);

    // cout<<"found ratio"<<endl;


    double draw = uniformDistribution(generator);

    if(ratio>draw)
    {
        positions = suggestion;
        olddrift(dimension,particle) = drift;
        oldwavesquared = newwavesquared;
        return 1;
    }
    else return 0;
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
    newwavesquared = wavefunctionSquaredJastrow(suggestion);
    return newwavesquared/oldwavesquared;
}

double vmc::findRatio(int d, int p)
{
    newwavesquared = wavefunctionSquared(suggestion);
    return newwavesquared/oldwavesquared;
}

double vmc::findRatioImportance(int d, int p)
{
    return proposalratio(d,p)*findRatio(d,p);
}

double vmc::findRatioImportanceJastrow(int d, int p)
{
    return proposalratio(d,p)*findRatioJastrow(d,p);
}

double vmc::wavefunctionSquared(mat pos)
{
    return exp(-alpha*omega*accu(pos%pos));
}
double vmc::wavefunctionSquaredJastrow(mat pos)
{
    double exponent = 0;
    double r;
    for(int i = 1; i < nParticles; i++)
    {
        for(int j = 0; j < i; j++)
        {
            r = rDifference(pos,i,j);
            exponent += r/(1+beta*r);
        }
    }
    exponent *= 2;
    return exp(exponent)*wavefunctionSquared(pos);

}

double vmc::proposalratio(int d, int p)
{
    return exp(0.5*(drift*drift - olddrift(d,p)*olddrift(d,p))*stepLength2
               + positions(d,p)*(drift + olddrift(d,p))
               - suggestion(d,p)*(drift + olddrift(d,p)));
}

/*-----------------------------------------------------------------------------------*/

double vmc::rDifference(mat pos, int p, int q)
{
    double xdiff = pos(0,p) - pos(0,q);
    double ydiff = pos(1,p) - pos(1,q);
    return sqrt(xdiff*xdiff + ydiff*ydiff);
}


/*------------------------------------------------------------------------------------*/
double vmc::localEnergy()
{
    return 0.5*(1-alpha*alpha)*omega*omega*accu(positions%positions) + 2.0*alpha*omega;
}

double vmc::localEnergyJastrow()
{
    double sum = 0;
    double r;
    double b;


    r = rDifference(positions,0,1);
    b = 1.0/(1+beta*r);

    sum += b*b*(-b*b - 1.0/r + 2*beta*b + alpha*omega*r);

    sum += localEnergy();
    return sum;
}

double vmc::localEnergyInteraction()
{
    double sum = 0;

    sum += 1.0/rDifference(positions,0,1);



    sum += localEnergy();
    return sum;
}

double vmc::localEnergyJastrowInteraction()
{
    double sum = 0;
    double r;
    double b;


    r = rDifference(positions,0,1);
    b = 1.0/(1+beta*r);

    sum += 1.0/r;
    sum += b*b*(-b*b - 1.0/r + 2*beta*b + alpha*omega*r);

    sum += localEnergy();
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
    return exp(-omega*alpha*accu(pos%pos)/2);
}

double vmc::wavefunctionJastrow(mat pos)
{
    double r = rDifference(pos,1,0);
    return exp( -omega*alpha*accu(pos%pos)/2 + r/(1+beta*r) );
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
    return sum - alpha*omega*pos(d,p);
}
/*-------------------------------------------------------------------------*/


void vmc::writeEnergies(int nCycles, const char* filename)
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


void vmc::steepestDescent(int nCycles, double gamma)
{
    updatePointers();
    if((!useJastrow)||(useNumDiff)||(!useInteraction))
    {
        cout<<"Not implemented for this configuration"<<endl;
    }
    else
    {
        //cout<<"hi"<<endl;
        double oldalpha = -100;
        double oldbeta = -100;

        for(int tsc = 0, iter = 0; tsc<5; iter++)
        {
            oldalpha = alpha;
            oldbeta = beta;

            //cout<<1<<endl;
            updateOld();
            //cout<<2<<endl;

            double Esum = 0;
            double Asum = 0;
            double Bsum = 0;
            double EAsum = 0;
            double EBsum = 0;

            //cout<<3<<endl;

            double deltaE,dEda,dEdb;


            //cout<<4<<endl;


            for(int i = 0; i < nCycles; i++)
            {
                //cout<<"a";
                metropolisMove();
                //cout<<"b";
                deltaE = localEnergyJastrowInteraction();
                //cout<<"c";
                dEda = alphaDeriv();
               // cout<<"d";
                dEdb = betaDeriv();
               // cout<<"e";
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
            //cout<<alpha<<", "<<beta<<endl;
            if((abs(oldalpha - alpha) < 1e-4)&&(abs(oldbeta - beta) < 1e-4))
                tsc++;
            else
            {    tsc = 0;
                if(iter>30)
                {
                    gamma /= 10;
                    iter = 0;
                    //cout<<"Max iter reached new gamma is "<<gamma<<endl;
                }
            }
        }
    }
    cout<<"Steepest Descent completed."<<endl;
    cout<<"New parameters are:"<<endl;
    cout<<"alpha = "<<alpha<<", beta = "<<beta<<endl;
}

double vmc::alphaDeriv()
{
    return -0.5*omega*accu(positions%positions);
}

double vmc::betaDeriv()
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
            sum -= s*s;
        }
    }

    return sum;

}

void vmc::printResults()
{
    cout<< "Acceptanceratio: "<<AcceptanceRatio <<endl;


    cout<<"E, meanvar" <<endl;
    cout<<energy<<","<<","<< meanvar <<endl;
}

void vmc::updateOld()
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

void vmc::updatePointers()
{
    // Make sure function pointers match current booleans
    // 16 cases
    /*------------------------------------------------------------------------------*/

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

}
