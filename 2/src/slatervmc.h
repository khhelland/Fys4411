#ifndef SLATERVMC_H
#define SLATERVMC_H

#include <armadillo>

class slatervmc
{
public:
    slatervmc();
    slatervmc(int nParticles, double, double w, double b);

    void run(int,int);

    void writeEnergies(int,const char*);

    void steepestDescent(int,double);

    void printResults();


    bool useJastrow = true;
    bool useInteraction = true;
    bool useImportanceSampling = true;


private:

    void distributeParticles();
    void setupMatrices();

    void updatePointers();
    void updateOld();
    void updateSlaters(int);

    int metropolisMove();

    void findSuggestionUniform(int,int);
    void findSuggestionImportanceSamplingwithJastrow(int,int);
    void findSuggestionImportanceSamplingnoJastrow(int,int);


    double findRatio(int,int);
    double findRatioJastrow(int,int);
    double findRatioImportance(int,int);
    double findRatioImportanceJastrow(int,int);


    double proposalratio(int,int);


    double localEnergy();
    double localEnergyJastrow();
    double localEnergyInteraction();
    double localEnergyJastrowInteraction();

    double kinetic();
    double kineticJastrow();

    double potential();
    double potentialInteraction();

    double slaterlaplace();
    double jastrowlaplace();

    arma::vec jastrowgrad(int);
    arma::vec slatergrad(int);


    double rDifference(arma::mat,int,int);

    double drifttermwithJastrow(arma::mat,int,int);
    double drifttermnoJastrow(arma::mat,int,int);

    double betaDeriv();

    std::default_random_engine generator;
    std::normal_distribution<double> normalDistribution;
    std::uniform_real_distribution<double> uniformDistribution;


    arma::mat positions;
    arma::mat suggestion;

    arma::mat slaterUp;
    arma::mat slaterDown;
    arma::mat inverseUp;
    arma::mat inverseDown;


    int nParticles;
    int nOrbitals;

    double stepLength;
    double stepLength2;


    double omega;
    double beta = 1;

    double a(int i,int j){return ((i<(nParticles/2))==(j<(nParticles/2))) ? 1.0/3.0 : 1.0;}

    double drift = 0;
    arma::mat olddrift;
    double oldwavesquared;
    double newwavesquared;

    //results
    double energy;
    double energySquared;
    double AcceptanceRatio;

    //function pointers;
    double (slatervmc::*localEnergyPointer)() = nullptr;
    void (slatervmc::*findSuggestionPointer)(int,int) = nullptr;
    double (slatervmc::*findRatioPointer)(int,int) = nullptr;
};

#endif // SLATERVMC_H
