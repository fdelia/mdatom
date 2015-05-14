/*-----------------FILE RANDOM.CPP --------------------------------*/

/* R. GEURTSEN, GRONINGEN, WFVG, JULY 1987
 * H.J.W. Spoelder, Amsterdadm, C-vserion 1.0, FEB 1988
 *
 *     RANDOM(RAND, SEED)
 *
 *     RANDOM GENERATES A RANDOM NUMBER RAND, USING A LINEAR
 *     CONGRUENTIAL METHOD. THE RECURSION FORMULA
 *
 *         IRAND = MOD(IRAND * B + 1, A)
 *
 *     IS USED WITH  B = 31415821  AND  A = 100000000. THE LAST DIGIT
 *     FROM THE RANDOM INTEGER IRAND IS CHOPPED OF, AND THE NUMBER
 *     IS SCALED TO A REAL VALUE RAND BETWEEN 0 AND 1, INCLUDING 0 BUT
 *     EXCLUDING 1.
 *
 *     RAND = DELIVERED WITH RANDOM NUMBER BETWEEN 0 AND 1
 *     IG = RANDOM NUMBER GENERATOR SEED, IS DELIVERED WITH RANDOM
 *          INTEGER
 *
 *
 */
/* 12-mar-91 had to add void declaration for a_random */
const int       M      =  100000000;
const int       M1     =  10000;
const int       MULT   =  31415821;

#define         abs(x)  (((x) > 0) ? (x) : -(x))


double a_random(int *seed){
    int             irand;
    int             irandh, irandl, multh, multl;
    double          r;

	irand = abs(*seed) % M;
    irandh = (int)(irand/M1);
    irandl = irand%M1;
    multh  = (int)(MULT/M1);
    multl  = MULT%M1;
    irand = ((irandh*multl + irandl*multh)%M1)*M1 + irandl*multl;
    irand = (irand+1)%M;
	r = (double)(irand/10.)*10./(double)M;
    if(r <= 0. || r > 1.) r = 0.;
    *seed = irand;
	return r;
}
