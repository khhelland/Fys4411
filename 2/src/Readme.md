Source files for project 2.
There are 2 .pro files

-2.pro is non-paralell. It includes
 -tests
 -test
 -main
 -vmc
 -slatervmc
-para2.pro is paralellized. It includes
 -paramain
 -slatervmc
 -vmc
-vmc is the two particle integrator class
-slatervmc is the N particle integrator class. It includes
 -ho2d, which  holds functions related to the single particle wavefunctions.
-blocking.py is a python script for blocking
-test.cpp holds some functions for benchmarks etc
