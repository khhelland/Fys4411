#ifndef VMC_H
#define VMC_H
#include <armadillo>
#include <random>


class vmc
{
public:
    vmc();
    vmc( double, double , double , double ,int );

    void run(int,int);

    void writeEnergies(int,const char*);

    void steepestDescent(int,double);

    void printResults();


    bool useJastrow = true;
    bool useInteraction = true;
    bool useNumDiff = false;
    bool useImportanceSampling = true;

    //results
    double energy;
    double energySquared;
    double error;
    double AcceptanceRatio;

    double mean_distance;


private:

    void distributeParticles();

    void updatePointers();
    void updateOld();

    int metropolisMove();

    void findSuggestionUniform(int,int);
    void findSuggestionImportanceSamplingwithJastrow(int,int);
    void findSuggestionImportanceSamplingnoJastrow(int,int);


    double findRatio(int,int);
    double findRatioJastrow(int,int);
    double findRatioImportance(int,int);
    double findRatioImportanceJastrow(int,int);

    double wavefunctionSquared(arma::mat pos);
    double wavefunctionSquaredJastrow(arma::mat pos);
    double proposalratio(int,int);


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

    double alphaDeriv();
    double betaDeriv();

    std::default_random_engine generator;
    std::normal_distribution<double> normalDistribution;
    std::uniform_real_distribution<double> uniformDistribution;


    arma::mat positions;
    arma::mat suggestion;

    int nParticles = 2;

    double stepLength;
    double stepLength2;

    double omega;
    double alpha;
    double beta = 1;

    double drift = 0;
    arma::mat olddrift;
    double oldwavesquared;
    double newwavesquared;



    //(void*)(int,int) function;
    double (vmc::*localEnergyPointer)() = nullptr;
    void (vmc::*findSuggestionPointer)(int,int) = nullptr;
    double (vmc::*findRatioPointer)(int,int) = nullptr;

};

#endif // VMC_H
