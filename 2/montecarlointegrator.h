#ifndef MONTECARLOINTEGRATOR_H
#define MONTECARLOINTEGRATOR_H


class montecarlointegrator
{
public:
    montecarlointegrator();
    montecarlointegrator(double,double,int);

    void move_walkers();
    double pdf(double x, double y);

    void integrate();

private:
    default_random_engine generator;
    normal_distribution<double> distribution;
    uniform_real_distribution<double> uniformDistribution;
    vector< vector< double> > walkers;

};

#endif // MONTECARLOINTEGRATOR_H
