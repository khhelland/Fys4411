#include <vector>
#include <random>

//sketch of metropolis program
using namespace std;


void move_particles(&vector < vector < double > > vec)
{
    default_random_engine generator;
    normal_distribution<double> distribution(0,1);
    uniform_real_distribution<double> uniformDistribution(0.0,1.0);
    double stepLength = 1;


    for(auto it = vec.begin(); it!= vec.end(); it++)
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
