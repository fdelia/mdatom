/*--------------------FILE FORCEA.C  #3 of 10------------------------------*/
/* 
 * w.f. van gunsteren, groningen, feb. 1988
 * h.j.w. spoelder, Amsterdam, C-version 1.0, feb. 1988
 * gee, 10/01, C++ version 1.0 .
 * 
 *     forcea (nat, x, box, epslj, siglj, rcutf, epot, f, vir,
 *             rcutg,ngr,igr)
 *
 *     forcea computes the potential energy , virial and
 *     interatomic forces for a lennard-jones potential
 *
 *         v(rij)  =  4*epslj*((siglj/rij)**12 - (siglj/rij)**6)
 *
 *     for all pairs (i<j) of atoms i and j at a distance rij.
 *     the pairs are selected using a cut-off radius rcutf.
 *     rectangular periodic boundary conditions are applied, that is,
 *     only nearest periodic images are included in the interaction.
 *     the contribution of the atom pairs to the radial distribution
 *     function g(r) (not yet properly scaled with r**2, etc.) is
 *     calculated too.
 *
 *     nat = number of atoms
 *     x(1..3*nat) = atom cartesian coordinates
 *     box(1..3) = lengths of the edges of the (periodic) box (>0)
 *     epslj = lennard-jones parameter epsilon
 *     siglj = lennard-jones parameter sigma
 *     rcutf = cut-off radius for selecting the atom pairs
 *             when calculating the interaction
 *     epot = delivered with the potential energy
 *     f(1..3*nat) = delivered with the forces
 *     vir = delivered with the virial
 *     rcutg = maximum distance rij for which the contribution to the
 *             radial distribution function g(rij) is evaluated
 *     ngr = number of intervals on the abcis of g(rij) between rij=0
 *           and rij=rcutg; one g(rij) value per interval is calculated
 *           (>0, rcutg>0)
 *     igr(1..ngr) : igr(n) is delivered with the number of ordered
 *                   atom pairs (i,j) at distance rij<rcutg, for which
 *                   n = rij/(rcutg/ngr) + 1, added to igr(n) (rcutg>0)
 *
 *
 */
#include        <cmath>    // For sqrt() .
#include <iostream>

using namespace std;

#define         nint(x)         (((x) > 0) ? (int)((x)+.5) : (int)((x)-.5))

void forcea( int nat, double x[], double box[], double epslj, double siglj,
             double rcutf, double *epot, double f[], double *vir,
             double rcutg, int ngr, int igr[]){

        double  boxinv[3];   // Inverse box length
        double  xij[3];      // Inter-particle vector
        double  rij;         // Inter-particle distance
        double  rij2;        // Inter-particle distance squared
        double  riji2;       // inverse inter-particle distance squared
        double  riji6;       // inverse inter-particle distance (6th power)

        double  sig6;        // L-J parameter
        double  c6;          // L-J parameter
        double  c12;         // L-J parameter
        double  crh;         // L-J parameter
        double  crhh;        // L-J parameter
        double  eij;         // L-J potential (work variable)

	double  rcutf2;      // long-range interaction cut-off (squared)
        double  rcutg2;      // radial distribution function cut-off (squared)

        double dij;          // Magnitude of force (derived from potential)
                             //  divided by inter-particle distance 
	double  df;          // Force increment in direction of inter-particle
                             //  vector 
                             //(note: xij[m]/rij is unit vector in
                             // inter-particle direction.)


        double  dgr;         // bin width for partitioning radial distance
                             //  in calculation of radial distribution 
                             //   function
        int     nat3;        // nr. degrees of freedom

        double one=1.0;
        int     m, i, j, i3, j3, n;  // Array counters
/*
 * initialise variables
 */
        nat3 = 3*nat;
        for(m=0; m<3; m++)
                boxinv[m] = one/box[m];
        sig6 = siglj*siglj;
        sig6 = sig6*sig6*sig6;
        c6   = 4.*epslj*sig6;
        c12  = c6*sig6;
        rcutf2 = rcutf*rcutf;
        rcutg2 = rcutg*rcutg;
        if(rcutg > 0)
                dgr = rcutg/ngr;
/*
 * set potential energy, virial and forces equal to zero
 */
        *epot = 0;
        *vir  = 0;
        for(j3=0; j3<nat3; j3++)
                f[j3] = 0;
/*
 * for all ordered pairs of atoms i and j
 */
        for(i=0; i<nat-1; i++) {
                i3 = 3*i;
                for(j=i+1; j<nat; j++) {
                        j3 = 3*j;
/*
 * calculate distance rij and apply periodic boundary conditions
 */
                        rij2 = 0.0;
                        for(m=0; m<3; m++) {
                                xij[m] = x[i3+m]-x[j3+m];
                                xij[m] = xij[m]-nint(xij[m]*boxinv[m])*box[m];
                                rij2  += xij[m]*xij[m];
                        }
/*
 * calculate potential energy, forces and virial
 */
                        if(rij2 < rcutf2) {
                                riji2  = one/rij2;
                                riji6  = riji2*riji2*riji2;
                                crh    = c12*riji6;
                                crhh   = crh-c6;
                                eij    = crhh*riji6;
                                *epot += eij;
                                dij   = 6.*(crh+crhh)*riji6*riji2;
                                for(m=0; m<3; m++) {
                                        df = xij[m]*dij;
                                        f[i3+m] += df;
                                        f[j3+m] -= df;
                                        *vir -= xij[m]*df;
                                }
                        }
/*
 * evaluate contribution to the radial distribution function
 */
                        if(rij2 < rcutg2) {
                                rij = sqrt(rij2);
                                n   = (int)(rij/dgr);
                                n   = (n > ngr) ? ngr : n;
                                igr[n]++;
                        }
                }
        }
        *vir /= 2.;
}






