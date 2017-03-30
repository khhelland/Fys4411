#ifndef VMC_H
#define VMC_H
#include <armadillo>
#include <random>


class vmc
{
public:
    vmc();
    vmc(int nParticles, double rSigma);



//private:

    void run(int);

    double wavefunctionSquared(arma::mat);
    double localEnergy(arma::mat);

    std::default_random_engine generator;
    std::normal_distribution<double> distribution;
    std::uniform_real_distribution<double> uniformDistribution;


    arma::mat positions;
    arma::mat suggestion;

    int nParticles;




};

#endif // VMC_H
