/*-------------------FILE SHIA.C   #10 of 10------------------------------*/

/*
 * W.F. VAN GUNSTEREN, GRONINGEN, FEB. 1988
 * H.J.W. Spoelder, Amsterdam, C-Version 1.0, FEB. 1988
 *
 *     SHIA(NAT, X, BOX, XMIN)
 *
 *     SHIA TRANSLATES ATOMS, APPLYING PERIODIC BOUNDARY CONDITIONS
 *     SUCH THAT THE ATOMS WILL LIE IN THE SPECIFIED (PERIODIC) BOX.
 *     THIS BOX LIES IN THE POSITIVE QUADRANT WITH RESPECT TO AN
 *     ORIGIN XMIN(1..3).
 *
 *     NAT = NUMBER OF ATOMS
 *     X(1..3*NAT) = ATOM CARTESIAN COORDINATES
 *     BOX(1..3) = LENGTHS OF THE EDGES OF THE (PERIODIC) BOX (>0)
 *     XMIN(1..3) = CARTESIAN COORDINATES OF THE ORIGIN OF THE BOX
 *
 *
 */
#define         nint(x)         (((x) > 0) ? (int)((x)+.5) : (int)((x)-.5))

void shia(int nat, double x[], double box[], double xmin[]){

        double  boxinv[3], xmb[3], one=1., xh;
        int     i, j, j3;

        for(i=0; i<3; i++) {
                boxinv[i] = one/box[i];
                xmb[i] = xmin[i]+box[i]/2.;
        }
        for(j=0; j<nat; j++) {
                j3 = 3*j;
                for(i=0; i<3; i++) {
                        xh = (x[j3+i]-xmb[i])*boxinv[i];
                        x[j3+i] -= nint(xh)*box[i];
                }
        }
}
