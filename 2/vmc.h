#ifndef VMC_H
#define VMC_H
#include <armadillo>
#include <random>


class vmc
{
public:
    vmc();
    vmc(int nParticles, double, double alph, double w, double b);

    void run(int,int);

    bool useJastrow = true;
    bool useInteraction = true;
    bool useNumDiff = false;
    bool useImportanceSampling = true;


//private:




    void distributeParticles();

    double wavefunctionSquared(arma::mat);
    double localEnergy(arma::mat);

    double rDifference(arma::mat,int,int);
    double laplacian(arma::mat);
    double wavefunction(arma::mat);

    double driftterm(arma::mat, int, int dim);
    double proposalDensity(double,double,double);

    std::default_random_engine generator;
    std::normal_distribution<double> normalDistribution;
    std::uniform_real_distribution<double> uniformDistribution;


    arma::mat positions;
    arma::mat suggestion;

    int nParticles;

    double stepLength;
    double stepLength2;

    double alpha;
    double omega;
    double beta = 1;
    double a = 1;


};

#endif // VMC_H
