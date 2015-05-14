/*------------------FILE CONFA.C  #2 of 10-----------------------------*/
/* w.f. van gunsteren, groningen, feb. 1988
 * h.j.w. spoelder, Amsterdam, C-version 1.0, feb. 1988
 *
 *     confa(box, nbox, xmin, dx, ig, x)
 *
 *     confa puts nat = nbox(1)*nbox(2)*nbox(3) atoms on a periodic
 *     lattice, such that along each side (m=1..3) of length box(m)
 *     will lie nbox(m) atoms.
 *     the box lies in the positive quadrant with respect to an
 *     origin xmin(1..3).
 *     uniformly distributed random displacements of maximum size dx/2
 *     can be added to the atomic positions in each direction.
 *
 *     box(1..3) = lengths of the edges of the (periodic) box (>0)
 *     nbox(1..3) = number of atoms to be placed along each box edge
 *                  (>0)
 *     xmin(1..3) = cartesian coordinates of the origin of the box
 *     dx = maximum spread of an atomic position around a lattice point
 *     ig = random number generator seed
 *     x(1..3*nat) = delivered with the atom cartesian coordinates
 *
 *     confa uses subr. random.
 *
 */
/*
 * who am I to stuff around in here ?
 * Added void declaration 12-mar-91
 * Added prototype for random. Surely this is dangerous given the
 * number of systems which define random in some other way, especially
 * since here it is a void function !
 */
#include        "mdprog.h"

void confa(double box[], int nbox[], double xmin[], double dx, int *ig, double x[]){

        double  xh[3], dr, rand;
        int     nh[3], m, j, nx, ny, nz;

        for(j=0; j<3; j++) {
                xh[j] = 0.0;
                nh[j] = 0;
        }
/*
 * initialise variables
 */
        for(m=0; m<3; m++)
                xh[m] = box[m]/nbox[m];
/*
 * put the atoms on a lattice and add random displacements
 */
        j = 0;
        for(nx=0; nx<nbox[0]; nx++) {
                nh[0] = (nx+1);
                for(ny=0; ny<nbox[1]; ny++) {
                        nh[1] = (ny+1);
                        for(nz=0; nz<nbox[2]; nz++) {
                                nh[2] = (nz+1);
                                for(m=0; m<3; m++) {
                                        rand = a_random(ig);
                                        dr = (rand-0.5)*dx;
                                        x[j+m] = xmin[m]+nh[m]*xh[m]+dr;
                                }
				j += 3;
                        }
                }
        }
}
