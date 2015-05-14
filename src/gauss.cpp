/*--------------------FILE GAUSS.C  #4 of 10-------------------------------*/
/* W.F. VAN GUNSTEREN, GRONINGEN, JULY 1987
 * H.J.W. Spoelder, Amsterdam, C-version 1.0, Feb 1988
 *
 *     GAUSS(AM, SD, V, SEED)
 *
 *     GAUSS WILL SUPPLY A NORMALLY DISTRIBUTED RANDOM NUMBER WITH
 *     A GIVEN MEAN AND STANDARD DEVIATION. IT USES 12 UNIFORM RANDOM
 *     NUMBERS TO COMPUTE NORMAL RANDOM NUMBERS BY THE CENTRAL LIMIT
 *     THEOREM (IBM 360 SCI.SUBR.PACKAGE).
 *
 *     AM   = THE DESIRED MEAN OF THE NORMAL DISTRIBUTION
 *     SD   = THE DESIRED STANDARD DEVIATION OF THE NORMAL DISTRIBUTION
 *     V    = VALUE OF THE COMPUTED NORMAL RANDOM VARIABLE
 *     SEED = RANDOM NUMBER GENERATOR SEED
 *
 *     GAUSS USES AS RANDOM NUMBER GENERATOR SUBR. RANDOM(Y,IG)
 *     (GROMOS, ANY MACHINE), OR SUBR. RANDU(IG,IY,Y) ON IBM MACHINES,
 *     FUNCTION RAN(IG) ON VAX MACHINES,
 *     OR INTRINSIC FUNCTION RANF(N) ON CDC MACHINES.
 *
 *
 */
#include "mdprog.h"
double gauss(double am, double sd, int *seed){

        double  a = 0;
        int     i;

        for(i=0; i<12; i++) {
                a += a_random(seed);
        }
        return (a-6.)*sd+am;
}
