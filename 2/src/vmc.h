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

    int blocksize(int);

    void printResults();


    bool useJastrow = true;
    bool useInteraction = true;
    bool useNumDiff = false;
    bool useImportanceSampling = true;


private:

    void distributeParticles();

    void findSuggestionUniform(int,int);
    void findSuggestionImportanceSamplingwithJastrow(int,int);
    void findSuggestionImportanceSamplingnoJastrow(int,int);


    double findRatio(int,int);
    double findRatioJastrow(int,int);
    double findRatioImportance(int,int);
    double findRatioImportanceJastrow(int,int);

    double localEnergy();
    double localEnergyJastrow();
    double localEnergyInteraction();
    double localEnergyNumdiff();
    double localEnergyJastrowInteraction();
    double localEnergyJastrowNumdiff();
    double localEnergyInteractionNumdiff();
    double localEnergyJastrowInteractionNumdiff();

    double wavefunction(arma::mat);
    double wavefunctionJastrow(arma::mat);

    double rDifference(arma::mat,int,int);

    double drifttermwithJastrow(arma::mat,int,int);
    double drifttermnoJastrow(arma::mat,int,int);


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

    double drift = 0;
    arma::mat olddrift;

    //results
    double energy;
    double energySquared;
    double AcceptanceRatio;



};

#endif // VMC_H