//#include "integrator.h"
//#include <harmonicoscillator2d.h>
//#include <math.h>
//#include <omp.h>
//#include <iostream>

////// The integrand in Cartesian coordinates
////double integrator::Integrand(double, double,double, double);
////// The Gauss Hermite integration function
////double  integrator::GaussHermiteIntegration(int);

////// Getting the Gaussian quadrature weights and integration points
////void integrator::GaussHermiteQuadrature(double *, double *, int);

//using namespace std;


//integrator::integrator(HarmonicOscillator2D a ,HarmonicOscillator2D b ,HarmonicOscillator2D c ,HarmonicOscillator2D d)
//{
//   p = a;
//   q = b;
//   r = c;
//   s = d;

//}


//double integrator::integrate(int n)
//{
//    int thread_num;
//    double wtime;
//    omp_set_num_threads(8);
//    thread_num = omp_get_max_threads ();
//    wtime = omp_get_wtime ( );
//    double Integral = GaussHermiteIntegration(n);
//    /*wtime = omp_get_wtime ( ) - wtime;
//    cout << setiosflags(ios::showpoint | ios::uppercase);
//    cout << setprecision(15) << setw(20) << "Time used  integration=" << wtime  << endl;
//    cout << "Final Integral    = " << Integral << endl;*/
//    return Integral;
//}



//double integrator::Integrand(double x1, double y1, double x2, double y2)
//{
//    if ((x1==x2) && (y1==y2)) { return 0;}
//    double coulomb = 1;//(sqrt(pow(x1-x2,2)+pow(y1-y2,2)));

//    return p.wavefunction_no_exp(x1,y1)*q.wavefunction_no_exp(x1,y1)*r.wavefunction_no_exp(x2,y2)*s.wavefunction_no_exp(x2,y2)*coulomb;
//}


////Plain Gauss-Hermite integration with cartesian variables, brute force
//double  integrator::GaussHermiteIntegration(int n)
//{
//  double *x = new double [n];
//  double *w = new double [n];
//  GaussHermiteQuadrature(x, w, n);
//  double Integral = 0.0;
//  int i, j, k, l;
//  # pragma omp parallel default(shared) private (i, j, k, l) reduction(+:Integral)
//  {
//    # pragma omp for
//    for (i = 0;  i < n; i++)
//    {
//        for (j = 0;  j < n; j++)
//        {
//            for (k = 0;  k < n; k++)
//            {
//                for (l = 0;  l < n; l++)
//                {
//                    Integral += w[i]*w[j]*w[k]*w[l]*Integrand(x[i],x[j],x[k],x[l]);
//                }
//            }
//        }
//    }
//  } // end parallel region
//  delete [] x;
//  delete [] w;
//  return Integral;
//}



//// Setting Gaussian quadrature weights and integration points
//void integrator::GaussHermiteQuadrature(double *x, double *w, int n)
//{
//  int i,its,j,m;
//  double p1,p2,p3,pp,z,z1;
//  double Epsilon = 3.0e-14, PIM4 = 0.7511255444649425;
//  int MaxIterations = 10;
//  m=(n+1)/2;
//  for (i=1;i<=m;i++) {
//    if (i == 1) {
//      z=sqrt((double)(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667);
//    } else if (i == 2) {
//      z -= 1.14*pow((double)n,0.426)/z;
//    } else if (i == 3) {
//      z=1.86*z-0.86*x[0];
//    } else if (i == 4) {
//      z=1.91*z-0.91*x[1];
//    } else {
//      z=2.0*z-x[i-3];
//    }
//    for (its=1;its<=MaxIterations;its++) {
//      p1=PIM4;
//      p2=0.0;
//      for (j=1;j<=n;j++) {
//    p3=p2;
//    p2=p1;
//    p1=z*sqrt(2.0/j)*p2-sqrt(((double)(j-1))/j)*p3;
//      }
//      pp=sqrt((double)2*n)*p2;
//      z1=z;
//      z=z1-p1/pp;
//      if (fabs(z-z1) <= Epsilon) break;
//    }
//    if (its > MaxIterations) cout << "too many iterations in Hermite quadrature" << endl;
//    x[i-1]=z;
//    x[n-i] = -z;
//    w[i-1]=2.0/(pp*pp);
//    w[n-i]=w[i-1];
//  }
//} // Gaussian quadrature weights and integration points

