#ifndef HO2D_H
#define HO2D_H
#include <armadillo>
double ho2d(int,double,double,double);
arma::vec ho2dgrad(int deg, double w, double x, double y);
double ho2ddw(int deg, double w, double x, double y);

arma::vec quantumnumbers(int);
double hermite(int,double);
double hermitederiv(int d, double x, int n);
double ho2denergy(int deg);

#endif // HO2D_H
