#include <harmonicoscillator2d.h>
#include <coulomb_functions.h>
#include <hartreefock.h>
#include <armadillo>


using namespace arma;



hartreefock::hartreefock(int particles , int shells , double w):
    shells(shells), n_particles(particles), frequency(w)
{
    n_orbitals = shells*(shells+1);
    find_ref_energies();
    hf_energies = zeros<vec>(n_orbitals);
    old_energies = zeros<vec>(n_orbitals);
    C = eye<mat>(n_orbitals,n_orbitals);
    densitymatrix = zeros<mat>(n_orbitals,n_orbitals);
    updatedensitymatrix();
    fockmatrix = zeros<mat>(n_orbitals,n_orbitals);

}

void hartreefock::find_ref_energies()
{
    ref_energies = zeros<vec>(n_orbitals);
    for(int i = 0; i < n_orbitals; i++)
    {
        HarmonicOscillator2D a(i);
        ref_energies(i) = a.getenergy();
    }
}


double hartreefock::matrixelement_as(int ip,int iq, int ir, int is)
{
    HarmonicOscillator2D p(ip);
    HarmonicOscillator2D q(iq);
    HarmonicOscillator2D r(ir);
    HarmonicOscillator2D s(is);

    //spin and m conservation
    if ((p.spin + q.spin != r.spin + s.spin)&&(p.m+q.m != r.m + s.m))
    {
        return 0;
    }


    double direct = Coulomb_HO(frequency,p.n,p.m,q.n,q.m,r.n,r.m,s.n,s.m);
    double exchange = Coulomb_HO(frequency,p.n,p.m,q.n,q.m,s.n,s.m,r.n,r.m);
    return direct-exchange;
}



void hartreefock::updatefockmatrix()
{
    for(int alpha = 0;alpha < n_orbitals;alpha++)
    {
        for(int beta = 0; beta <= alpha; beta++)
        {
            double s = 0;
            for(int gamma = 0;gamma < n_orbitals;gamma++)
            {
                for(int delta = 0; delta < n_orbitals; delta++)
                {
                    s+= densitymatrix(gamma,delta)*matrixelement_as(alpha,gamma,beta,delta);
                }
             }
            if (alpha == beta){s+=ref_energies(alpha);}
            fockmatrix(alpha,beta) = s;
            fockmatrix(beta,alpha) = s;
        }
    }
}

void hartreefock::updatedensitymatrix()
{
    int s;
    for(int gamma = 0;gamma < n_orbitals;gamma++)
    {
        for(int delta = 0; delta <= gamma; delta++)
        {
            s = 0;
            for(int i = 0; i < n_particles; i++)
            {
                s += C(i,gamma)*C(i,delta);
            }
            densitymatrix(gamma,delta) = s;
            densitymatrix(delta,gamma) = s;
        }
    }
}


void hartreefock::diagonalizefockmatrix()
{
    eig_sym(hf_energies,C,fockmatrix);
}


double hartreefock::finddifference()
{
    return sum(abs(hf_energies-old_energies))/n_orbitals;
}


void hartreefock::run(int maxcount, double epsilon)
{   double difference = 1;
    for(int count = 0;(count<=maxcount)&&(difference>epsilon);count++)
    {
        old_energies = hf_energies;
        updatedensitymatrix();
        updatefockmatrix();
        diagonalizefockmatrix();
        difference = finddifference();
    }
}

double hartreefock::getenergy()
{
    return sum(hf_energies);
}
