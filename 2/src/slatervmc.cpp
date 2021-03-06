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
 *Class for many particle MC integration.
 * Use: initialize, change booleans from default if needed
 * and run.
 *
 */


slatervmc::slatervmc()
{

}
slatervmc::slatervmc(int nParticles, double step, double w, double a, double b, int s):
    nParticles(nParticles),
    stepLength(step),
    omega(w),
    alpha(a),
    beta(b)
{
    normalDistribution = normal_distribution<double>(0.0,1.0);
    uniformDistribution = uniform_real_distribution<double>(0.0,1.0);
    nOrbitals = nParticles/2;
    alphaomega = alpha*omega;
    distributeParticles();
    setupMatrices();
    olddrift = mat(2,nParticles);
    stepLength2 = stepLength*stepLength;
    generator.seed(s);

}

void slatervmc::distributeParticles()
{
    //chance of coincident particles should be ~ 0
    positions = mat(2,nParticles,fill::randu);

//    for(int p = 0;p<nParticles;p++)
//    {
//        positions(0,p) = p;
//        positions(1,p) = p;
//    }
}

void slatervmc::setupMatrices()
{
    slaterDown = mat(nOrbitals,nOrbitals);
    slaterUp = mat(nOrbitals,nOrbitals);

    for(int i = 0; i < nOrbitals; i++)
    {
        for(int j = 0; j<nOrbitals; j++)
        {
            slaterDown(i,j) = ho2d(j,alphaomega,positions(0,i),positions(1,i));
            slaterUp(i,j) = ho2d(j,alphaomega,positions(0,i+nOrbitals),positions(1,i+nOrbitals));
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
    for(int i = 0; i<2*blocksize ; i++) metropolisMove();

    double Esum = 0;
    double Esum2 = 0;
    double Eblock = 0;

    double Ksum = 0;
    double Ksum2 = 0;
    double Kblock = 0;
    double Vsum = 0;
    double Vsum2 = 0;
    double Vblock = 0;


    int acceptcount = 0;
    for(int i = 1;i<=nCycles;i++)
    {
        acceptcount += metropolisMove();

        Kblock += (this->*kineticPtr)();
        Vblock += (this->*potentialPtr)();


        if(i % blocksize == 0)
        {
            Kblock /= blocksize;
            Vblock /= blocksize;
            Eblock = Kblock + Vblock;

            Ksum += Kblock;
            Ksum2 += Kblock*Kblock;
            Vsum += Vblock;
            Vsum2 += Vblock*Vblock;
            Esum += Eblock;
            Esum2 += Eblock*Eblock;

            //Eblock = 0;
            Kblock = 0;
            Vblock = 0;
        }
    }
    double dBlocks = (double)nCycles/blocksize;

    Esum /= dBlocks;
    Esum2 /= dBlocks;
    Ksum /= dBlocks;
    Ksum2 /= dBlocks;
    Vsum /= dBlocks;
    Vsum2 /= dBlocks;

    kinetic_energy = Ksum;
    kineticMeanVar = (Ksum2 - Ksum*Ksum)/(dBlocks-1);
    potential_energy = Vsum;
    potentialMeanVar = (Vsum2 - Vsum*Vsum)/(dBlocks -1);
    energy = Esum;
    energyMeanVar = (Esum2 -Esum*Esum)/(dBlocks - 1);


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
            R += ho2d(j, alphaomega, xnew, ynew)*inverseDown(j,p);
        }
    }
    else
    {
        for(int j = 0; j <nOrbitals; j++)
        {
            R += ho2d(j, alphaomega, xnew, ynew)*inverseUp(j,p-nOrbitals);
        }
    }
    return R*R;
}

double slatervmc::findRatioJastrow(int d ,int p)
{
    //funker for 2
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


double slatervmc::proposalratio(int d, int p)
{
    return exp(0.5*(drift*drift - olddrift(d,p)*olddrift(d,p))*stepLength2
               + positions(d,p)*(drift + olddrift(d,p))
               - suggestion(d,p)*(drift + olddrift(d,p)));
//    return exp(
//                (pow(suggestion(d,p) - positions(d,p) - stepLength2*olddrift(d,p),2)
//                - pow(positions(d,p) - suggestion(d,p) - stepLength2*drift,2)
//                 )/(2*stepLength2)
//                );
}

/*-----------------------------------------------------------------------------------*/

double slatervmc::rDifference(mat pos, int p, int q)
{
    double xdiff = pos(0,p) - pos(0,q);
    double ydiff = pos(1,p) - pos(1,q);
    return sqrt(xdiff*xdiff + ydiff*ydiff);
}



double slatervmc::kinetic()
{
    return -0.5*slaterlaplace();
}
double slatervmc::kineticJastrow()
{
    return -0.5*(slaterlaplace() + jastrowlaplace());
}
double slatervmc::potential()
{
    return 0.5*omega*omega*accu(positions%positions);
}
double slatervmc::potentialInteraction()
{
    double sum = 0;
    for(int i = 1; i < nParticles; i++)
    {
        for(int j = 0; j < i; j++)
        {
            sum += 1/rDifference(positions,i,j);
        }
    }
    return sum + potential();
}

double slatervmc::slaterlaplace()
{
    return alphaomega*alphaomega*accu(positions%positions) - 4*alphaomega*ho2denergy(nOrbitals);
}

double slatervmc::jastrowlaplace()
{
    double sum = 0;

    double av,b,r;
    for(int i = 1; i < nParticles ; i++)
    {
        for(int j = 0; j< i; j++)
        {
            r = rDifference(positions,i,j);
            b = 1/( 1 + beta*r);
            av = a(i,j);
            sum += av*b*b*(1/r - 2*beta*b);
        }
    }

    sum *= 2;

    for(int i = 0; i < nParticles; i++)
    {
        vec sgrad,jgrad;
        jgrad = jastrowgrad(i);
        sgrad = slatergrad(i);
        sum += dot(jgrad,jgrad) + 2*dot(jgrad,sgrad);

    }

    return sum;
}

vec slatervmc::slatergrad(int p)
{
    vec sum = zeros<vec>(2);
    if(p<nOrbitals)
    {
        for(int i = 0; i< nOrbitals; i++)
        {
            sum += ho2dgrad(i,alphaomega,positions(0,p),positions(1,p))*inverseDown(i,p);
        }
    }
    else
    {
        for(int i = 0; i< nOrbitals; i++)
        {
            sum += ho2dgrad(i,alphaomega,positions(0,p),positions(1,p))*inverseUp(i,p-nOrbitals);
        }
    }
    return  sum;
}


vec slatervmc::jastrowgrad( int p)
{
    //(1/J) dJ/dz(d)_p
    double sumx = 0;
    double sumy = 0;
    double av,r,b;
    for(int j = 0; j < nParticles; j++)
    {
        if(j != p)
        {
            av = a(p,j);
            r = rDifference(positions,j,p);
            b = 1/(1 + beta*r);
            av = av*b*b/r;
            sumx += av*(positions(0,p) - positions(0,j));
            sumy += av*(positions(1,p) - positions(1,j));
        }
    }
    vec res = {sumx,sumy};
    return res;
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
double slatervmc::drifttermnoJastrow(mat pos, int d , int p)
{
    return slatergrad(p)(d);
}

double slatervmc::drifttermwithJastrow(mat pos,int d, int p)
{
    return slatergrad(p)(d) + jastrowgrad(p)(d);
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
        e = (this->*kineticPtr)() + (this->*potentialPtr)();
        out<<e<<endl;
    }
    out.close();
}

void slatervmc::steepestDescentAuto(int nCycles, double gamma)
{
    updatePointers();
    vec params;
    if(useJastrow)
    {
        params = {alpha,beta};
    }
    else
    {
        params = {alpha};
    }
    vec grad;


    for(int tsc = 0, iter = 0; tsc<5; iter++)
    {

        grad = paramgrad(nCycles);

        while(alpha < gamma*grad(0))
        {
            gamma /= 10;
            cout<<"gamma too big, trying again with new gamma "<<gamma<<endl;
        }

        params -= gamma*grad;


        if(gamma*norm(grad) < 1e-5)
        {
            tsc++;
        }
        else
        {
            tsc = 0;

            if(iter>30)
            {
                gamma /= 10;
                iter = 0;
                cout<<"Max iter reached new  gamma is "<<gamma<<endl;
            }
        }

        updateParams(params);

        trans(params).print();

    }

    cout<<"Steepest Descent completed."<<endl;
    cout<<"New parameters are:"<<endl;
    cout<<"alpha = "<<alpha<<", beta = "<<beta<<endl;
}

void slatervmc::steepestDescentMan(int nCycles, double gamma)
{
    updatePointers();
    vec params;
    if(useJastrow)
    {
        params = {alpha,beta};
    }
    else
    {
        params = {alpha};
    }
    vec grad;

    cout<<"This function will run 1000 times if not stopped."<<endl;
    for(int i= 0;i<1000;i++)
    {
        trans(params).print();

        grad = paramgrad(nCycles);
        //trans(gamma*grad).print();
        params -= gamma*grad;
        updateParams(params);



    }
}



void slatervmc::updateParams(vec params)
{
    if(useJastrow)
    {
        beta = params(1);
    }
    alpha = params(0);
    alphaomega = alpha*omega;
    setupMatrices();
}

vec slatervmc::paramgrad(int nCycles)
{

    updateOld();

    double Esum = 0;
    double Asum = 0;
    double EAsum = 0;

    double deltaE,dEda,Bsum,EBsum,dEdb;

    if(useJastrow)
    {
        Bsum = 0;
        EBsum = 0;
    }



    for(int i = 0; i < nCycles; i++) metropolisMove();

    if(useJastrow)
    {
        for(int i = 0; i < nCycles; i++)
        {
            metropolisMove();
            deltaE = (this->*kineticPtr)() + (this->*potentialPtr)();

            dEda = alphaDeriv();
            dEdb = betaDeriv();

            Esum += deltaE;
            Asum += dEda;
            Bsum += dEdb;
            EAsum += deltaE*dEda;
            EBsum += deltaE*dEdb;

        }
    }
    else
    {
        for(int i = 0; i < nCycles; i++)
        {
            metropolisMove();
            deltaE = (this->*kineticPtr)() + (this->*potentialPtr)();

            dEda = alphaDeriv();

            Esum += deltaE;
            Asum += dEda;
            EAsum += deltaE*dEda;
        }
    }
    Esum /= nCycles;
    Asum /= nCycles;
    EAsum /= nCycles;

    vec grad;
    if(useJastrow)
    {
        Bsum /= nCycles;
        EBsum /= nCycles;

        grad = {2*(EAsum -Esum*Asum),2*(EBsum -Esum*Bsum)};
    }
    else
    {
        grad = {2*(EAsum - Esum*Asum)};
    }
    return grad;

}


double slatervmc::alphaDeriv()
{
    double sum = 0;
    for(int i = 0; i < nOrbitals; i++)
    {
        for(int j = 0; j < nOrbitals; j++)
        {
            sum +=ho2ddw(j,alphaomega,positions(0,i),positions(1,i))*inverseDown(j,i);
            sum +=ho2ddw(j,alphaomega,positions(0,i+nOrbitals),
                        positions(1,i+nOrbitals))*inverseUp(j,i);
        }
    }
    sum *= omega;
    //cout<<sum<<endl;
    return sum;
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


    cout<<"K\tsK\tV\tsV\tE\t sE" <<endl;
    cout<<kinetic_energy<<"\t"<<kineticMeanVar<<"\t"
        <<potential_energy<<"\t"<<potentialMeanVar<<"\t"
        <<energy<<"\t"<<energyMeanVar <<endl;
}

//double * slatervmc::getResults()
//{
//    double result[2];
//    result[0] = energy;
//    result[1] = energyError;
//    return result;
//}

void slatervmc::updateSlaters(int p)
{
    double xnew = positions(0,p);
    double ynew = positions(1,p);

    if(p<nOrbitals)
    {

        for(int j = 0; j<nOrbitals; j++)
        {
            slaterDown(p,j) = ho2d(j,alphaomega, xnew, ynew);
        }
        //Consider inserting efficient algorithm here
        inverseDown = inv(slaterDown);
    }
    else
    {
        for(int j = 0; j<nOrbitals; j++)
        {
            slaterUp(p-nOrbitals,j) = ho2d(j,alphaomega, xnew, ynew);
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
        kineticPtr = &slatervmc::kineticJastrow;
        potentialPtr = &slatervmc::potentialInteraction;
    }

    else if ((useImportanceSampling)&&(useInteraction)&&(!useJastrow) )
    {
        findSuggestionPointer = &slatervmc::findSuggestionImportanceSamplingnoJastrow;
        findRatioPointer = &slatervmc::findRatioImportance;
        kineticPtr = &slatervmc::kinetic;
        potentialPtr = &slatervmc::potentialInteraction;
    }

    else if ((useImportanceSampling)&&(!useInteraction)&&(useJastrow) )
    {
        findSuggestionPointer = &slatervmc::findSuggestionImportanceSamplingwithJastrow;
        findRatioPointer = &slatervmc::findRatioImportanceJastrow;
        kineticPtr = &slatervmc::kineticJastrow;
        potentialPtr = &slatervmc::potential;
    }

    else if ((!useImportanceSampling)&&(useInteraction)&&(useJastrow) )
    {
        findSuggestionPointer = &slatervmc::findSuggestionUniform;
        findRatioPointer = &slatervmc::findRatioJastrow;
        kineticPtr = &slatervmc::kineticJastrow;
        potentialPtr = &slatervmc::potentialInteraction;
    }


    else if ((!useImportanceSampling)&&(!useInteraction)&&(useJastrow) )
    {
        findSuggestionPointer = &slatervmc::findSuggestionUniform;
        findRatioPointer = &slatervmc::findRatioJastrow;
        kineticPtr = &slatervmc::kineticJastrow;
        potentialPtr = &slatervmc::potential;
    }

    else if ((!useImportanceSampling)&&(useInteraction)&&(!useJastrow) )
    {
        findSuggestionPointer = &slatervmc::findSuggestionUniform;
        findRatioPointer = &slatervmc::findRatio;
        kineticPtr = &slatervmc::kinetic;
        potentialPtr = &slatervmc::potentialInteraction;
    }

    else if ((useImportanceSampling)&&(!useInteraction)&&(!useJastrow) )
    {
        findSuggestionPointer = &slatervmc::findSuggestionImportanceSamplingnoJastrow;
        findRatioPointer = &slatervmc::findRatioImportance;
        kineticPtr = &slatervmc::kinetic;
        potentialPtr = &slatervmc::potential;
    }

    else if ((!useImportanceSampling)&&(!useInteraction)&&(!useJastrow) )
    {
        findSuggestionPointer = &slatervmc::findSuggestionUniform;
        findRatioPointer = &slatervmc::findRatio;
        kineticPtr = &slatervmc::kinetic;
        potentialPtr = &slatervmc::potential;
    }

    /*-------------------------------------------------------------------------------*/

}
