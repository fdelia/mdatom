/* 12-mar-91
 * Some of Wilfred's integers should have been parameter statements
 * in the original fortran.
 * They were converted into variables in the original C translation,
 * they are gathered here into #define's
 * posix compliant time functions 11/99 hschafer
 */
#ifndef MDPROG_H
#define MDPROG_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>


using namespace std;

const int MAXXCO = 10000;    /* number of degrees of freedom */
const int MAXNGR = 500;    /* maximum number of abcissa values for the  
                            * radial distribution 
                            */
const int MAXTITLE = 200;  /* Maximum number of characters for file titles */

const double one   = 1.0;

const string finalcoords_filename = "coords.final";
const string trajectory_filename = "coords.traj";

// Center of mass computation: cenma.cpp
void cenma( int nat, double x[], double v[], double amas, double *ekcm,
            double xcm[], double vcm[], int icm);

// Generation of initial atomic configuration: confa.cpp
void confa( double box[], int nbox[], double xmin[], double dx,
            int *ig, double x[]);

// Force, potential eneryg, virial and radial distribution function contributions: forcea.cpp
void forcea( int nat, double x[], double box[], double epslj, double siglj,
             double rcutf, double *epot, double f[], double *vir,
             double rcutg, int ngr, int igr[]);

// langevin.cpp
void langevin( int nat, double x[], double box[], double epslj, double siglj,
             double rcutf, double *epot, double f[], double *vir,
             double rcutg, int ngr, int igr[], double gamma, double amas, 
			 double boltz, double temp0, int ig, double v[]); // xxx

// Draw initial velocities from a Maxwellian distribution: gauss.cpp
double gauss ( double am, double sd, int *seed);

// Write time frames to disk: packa.cpp
void packa ( double r[], int ntot, int mode, string filename, int ntpw);

//void plot ( int nint, double gr[]);

// Obtain uniformly distributed random number: random.cpp
double a_random ( int *seed );

// Perform MD simulation: runmda.cpp
void runmda ( int nat, double x[], double v[], double f[], double amas,
              double epslj, double siglj, double rcutf, double box[],
              int nstlim, double t, double dt, int ntt, double temp0,
              double taut, double boltz, double *temp, int ntpr, int ntwx,
              int ntwxm, int ntpw, double rcutg, int ngr, int igr[] , double dtcoll, int ig, double gamma, double taup, double betat); // xxx

// Keep atoms in central periodic box: shia.cpp
void shia ( int nat, double x[], double box[], double xmin[] );

// Auxiliary functions: Implemented in auxiliary.cpp
void space (int n);
void usage  (void);
void readdbl(istream& in, double value[], int length);
void writedbl(ostream& out, double value[], int length);

#endif
