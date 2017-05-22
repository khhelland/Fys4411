#ifndef HO2D_H
#define HO2D_H
#include <armadillo>
double ho2d(int,double,double,double);
arma::vec ho2dgrad(int deg, double w, double x, double y);
double hermite(int,double);
double hermitederiv(int d, double x);
int ho2denergy(int nOrbitals);
#endif // HO2D_H
