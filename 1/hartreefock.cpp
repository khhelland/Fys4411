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
    fockmatrix = zeros<mat>(n_orbitals,n_orbitals);

    elementarray = new double[(int)pow(n_orbitals,4)];
}
hartreefock::~hartreefock()
{
    delete [] elementarray;
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

    double direct = 0;

    if((p.spin == r.spin)&&(q.spin == s.spin))
    {
        direct = Coulomb_HO(frequency,p.n,p.m,q.n,q.m,r.n,r.m,s.n,s.m);
    }

    double exchange = 0;

    if((p.spin == s.spin)&&(q.spin == r.spin))
    {
        exchange = Coulomb_HO(frequency,p.n,p.m,q.n,q.m,s.n,s.m,r.n,r.m);
    }
    return direct-exchange;
}



void hartreefock::updatefockmatrix()
{

    fockmatrix.zeros();

    for(int alpha = 0;alpha < n_orbitals;alpha++)
    {

        for(int beta = 0; beta <=alpha ; beta++)
        {
            double s = 0;
            for(int gamma = 0;gamma < n_orbitals;gamma++)
            {
                for(int delta = 0; delta < n_orbitals; delta++)
                {
                    s+= densitymatrix(gamma,delta)*elementarray[arrayindex(alpha,gamma,beta,delta)];
                            //matrixelement_as(alpha,gamma,beta,delta);
                }
             }
             fockmatrix(alpha,beta) = s;
             fockmatrix(beta,alpha) = s;
        }

        fockmatrix(alpha,alpha) += ref_energies(alpha);
    }
}

void hartreefock::updatedensitymatrix()
{
    //C.print();
    double s;
    for(int gamma = 0;gamma < n_orbitals;gamma++)
    {
        for(int delta = 0; delta <=gamma; delta++)
        {
            s = 0;
            densitymatrix(delta,gamma) = 0;
            for(int i = 0; i < n_particles; i++)
            {
                s += C(gamma,i)*C(delta,i);

            }
            densitymatrix(delta,gamma) = s;
            densitymatrix(gamma,delta) = s;
        }
    }


}


void hartreefock::diagonalizefockmatrix()
{
    eig_sym(hf_energies,C,fockmatrix);
    //C.print();
}


double hartreefock::finddifference()
{
    return sum(abs(hf_energies-old_energies))/n_orbitals;
}


double hartreefock::getenergy()
{
    double s = 0;
//    for(int i=0;i<n_particles; i++)
//    {
//        for(int j=0; j< n_particles; j++)
//        {
          for(int alpha = 0; alpha < n_orbitals; alpha++)
            {
                for(int beta = 0; beta< n_orbitals; beta++)
                {
                    for(int gamma = 0; gamma < n_orbitals; gamma++)
                    {
                        for(int delta = 0; delta< n_orbitals; delta++)
                        {
                            //s-= C(i,alpha)*C(i,beta)*C(j,gamma)*C(j,delta)*matrixelement_as(alpha,beta,gamma,delta);
                            s-= densitymatrix(alpha,gamma)*densitymatrix(beta,delta)
                                    *elementarray[arrayindex(alpha,beta,gamma,delta)];
                                    //*matrixelement_as(alpha,beta,delta,gamma);

                        }
                     }
                 }
             }
//         }
//    }

    s*=0.5;
    for(int i = 0; i< n_particles; i++)
    {
        s += hf_energies(i);
    }
    return s;
}

void hartreefock::print_sp_energy()
{
    hf_energies.print();
}


void hartreefock::findelements()
{
    for(int p = 0; p < n_orbitals; p++)
    {
        for(int q = 0; q < n_orbitals; q++)
        {
            for(int r = 0; r < n_orbitals; r++)
            {
                for(int s = 0; s < n_orbitals; s++)
                {
                    elementarray[arrayindex(p,q,r,s)] = matrixelement_as(p,q,r,s);
                    //cout<< elementarray[arrayindex(p,q,r,s)]<<" "<<matrixelement_as(p,q,r,s)<<endl;
                }
            }
        }
    }
}

int hartreefock::arrayindex(int p, int q , int r, int s)
{
    return p + q*n_orbitals + r*pow(n_orbitals,2) + s*pow(n_orbitals,3);
}



void hartreefock::run(int maxcount, double epsilon)
{
    findelements();
    double difference = 1;
    int count = 0;
    while((count<maxcount)&&(difference>epsilon))
    {
        old_energies = hf_energies;
        updatedensitymatrix();
        updatefockmatrix();
        diagonalizefockmatrix();
        difference = finddifference();
        count++;
    }
    if(count == maxcount)
    {
        cout<<"Warning: maximum number of iterations reached"<<endl;
    }
}
