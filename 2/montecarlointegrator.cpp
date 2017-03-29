#include "montecarlointegrator.h"
#include <math.h>


montecarlointegrator::montecarlointegrator()
{
    distribution(0,1);
    uniformDistribution(0.0,1.0);
}
montecarlointegrator::montecarlointegrator(double mu, double sigma, int N)
{
    distribution(mu,sigma);
    uniformDistribution(0.0,1.0);
    //create N 2d walkers at mu,mu
}

void montecarlointegrator::move_walkers()
{
    for(auto it = walkers.begin(); it!= walkers.end(); it++)
    {
        bool moved = 0;
        while(!moved)
        {
            double xstep = stepLength*distribution(generator);
            double ystep = stepLength*distribution(generator);
            ratio = pdf((*it)[0]+xstep,(*it)[1]+ystep)/pdf((*it)[0],(*it)[1]);
            if(ratio>uniformDistribution(generator))
            {
                (*it)[0] += xstep;
                (*it)[1] += ystep;
                moved = 1;
            }
        }
    }
}

double montecarlointegrator::pdf(double x,double y)
{
    return exp(-(pow(x,2)+pow(y,2)));
}
void montecarlointegrator::integrate()
{
    for(int i = 0;i<nSteps;i++)
        move_walkers();
}
