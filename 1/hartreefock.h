#ifndef HARTREEFOCK_H
#define HARTREEFOCK_H
#include <harmonicoscillator2d.h>
#include <armadillo>


class hartreefock
{
public:
    hartreefock(){};
    hartreefock(int particles ,int shells, double w);

    void run(int maxcount, double epsilon);
    double getenergy();


//private:
    int shells;
    int n_orbitals;
    int n_particles;
    double frequency;

    arma::mat densitymatrix;
    arma::mat C;
    arma::mat fockmatrix;

    arma::vec ref_energies;
    arma::vec hf_energies;
    arma::vec old_energies;

    void updatefockmatrix();
    void updatedensitymatrix();
    void diagonalizefockmatrix();
    void find_ref_energies();

    double matrixelement_as(int ip,int iq, int ir, int is);

    double finddifference();

private:

};




#endif // HARTREEFOCK_H

