/*---------------------FILE RUNMDA.CPP ----------------------------*/
/* w.f. van gunsteren, groningen, feb. 1988
 * h.j.w. spoelder, amsterdam, c-version 1.0, feb 1988
 * gee, 10/01, C++ version 1.0 .
 * mph, 01/15, removed tapes
 *
 *     function runmda (nat,x,v,f,amas,epslj,siglj,rcutf,box,
 *                      nstlim,t,dt,ntt,temp0,taut,boltz,temp,
 *                      ntpr,ntwx,ntwxm,ntpw,rcutg,ngr,igr)
 *
 *     runmda will carry out a molecular dynamics run for a
 *     collection of atoms.
 *     rectangular periodic boundary conditions are applied.
 *     it starts at time t, runs with time step dt and stops when the
 *     number of steps reaches nstlim.
 *     if ntt>0, the temperature t of the system is coupled to that of
 *     a bath of temperature temp0 with coupling time constant taut,
 *     see berendsen et al., j.chem.phys. 81 (1984) 3684.
 *     the differential equations are integrated using a leap-frog
 *     scheme.
 *
 *     function forcea supplies energies, forces and the virial.
 *     function shia keeps the atoms in the central periodic box.
 *     function cenma calculates the center of mass coordinates and motion.
 *     function packa is used for writing the atomic trajectories to disk.
 *
 *     results are printed every ntpr md steps.
 *     the center of mass motion is calculated every nlsq md steps.
 *     nlsq = max(nstlim/10,10,ntpr).
 *     finally, average results and fluctuations are printed.
 *
 *     nat = number of atoms
 *     x(1..3*nat) = initial atomic cartesian coordinates (at time t)
 *     v(1..3*nat) = initial atomic velocities at time t-dt/2
 *     x and v are delivered with final values (at time t and t-dt/2)
 *     f(1..3*nat) = work array (forces, etc.)
 *     amas = atomic mass (<>0)
 *     epslj = lennard-jones interaction parameter epsilon
 *     siglj = lennard-jones interaction parameter sigma
 *     rcutf = cut-off radius when calculating the interaction
 *     box(1..3) = lengths of the edges of the periodic box (<>0)
 *     nstlim = number of md-steps to be performed
 *     t = time at calling, delivered with final time
 *     dt = time step
 *     ntt = 0 : classical (constant total energy) md
 *         = 1 : md with velocity scaling (constant temperature)
 *     temp0 = reference temperature at which the system is kept (ntt>0)
 *     taut = temperature relaxation time (>0, ntt>0)
 *     boltz = proportionality constant (<>0) for calculating the
 *             temperature temp from the total kinetic energy ekin of
 *             the system using the relation: ekin = 3*nat*boltz*temp/2
 *     temp = delivered with average temperature
 *     ntpr : results (array ener(1..nren)) are printed every ntpr
 *            md steps
 *     ntwx,ntwxm : controls writing to tape12, see below
 *     ntpw : controls how function packa writes the data, see function packa
 *     rcutg = maximum distance r for which the radial distribution
 *             function g(r) is evaluated
 *     ngr = number of intervals on the abcis of g(r) between r=0 and
 *           r=rcutg; one g(r) value per interval is calculated
 *           (>0, rcutg>0)
 *     igr(1..ngr) : igr(n) is delivered with the total number of
 *                  ordered atom pairs at a distance r<rcutg for which
 *                   n = r/(rcutg/ngr)+1  (rcutg>0)
 *
 *     results can be written to trajectory_filename.
 *     atomic coordinates are written every ntwx steps
 *              for a maximum of ntwxm steps, using function packa.
 *
 *     runmda uses functions forcea, shia, cenma and packa.
 *
 */
#include        <iostream>    // For terminal I/O
#include        <iomanip>     // For terminal I/O format
#include        <cmath>       // For sqrt() .

#include        "mdprog.h"    // For run-time constants and 
                              //  function prototypes. 

using namespace std;

void runmda(int nat, double x[], double v[], double f[], double amas,
              double epslj, double siglj, double rcutf, double box[],
              int nstlim, double t, double dt, int ntt, double temp0,
              double taut, double boltz, double *temp, int ntpr, int ntwx,
              int ntwxm, int ntpw, double rcutg, int ngr, int igr[] , double dtcoll, int ig, double gamma, double taup, double betat){ // xxx

        double  xmin[3], xcm[3], vcm[3], ener[6], enert[6], enert2[6];
        double  one, fac, dtt, dt5, vh, pres, ekin0, vol, dtm, eold,
                ekg, epot, vir, scal, enew, vn, tspan, temp2, ekcm, lastcoll; // xxx
        int     nat3, i, j3, m, k, nhpr, nren, nstep, nlsq, collCounter; // xxx

/*; initialise */
        nat3 = 3*nat;
        one  = 1.;
        for(m=0; m<3; m++) {
	        xmin[m] = 0;
	    }
        fac = nat3*boltz/2.;
        ekin0 = fac*temp0;
        if(ntt > 0) {
	        dtt = dt/taut;
	    }
        dt5 = dt/2;
        dtm = dt/amas;
        vol = box[0]*box[1]*box[2];

        nhpr = 100*ntpr;
        nlsq = nstlim/10;
        if(nlsq < 10) {
	        nlsq = 10;
	    }
        if(nlsq < ntpr) {
	        nlsq = ntpr;
        }
	    if(nhpr > nlsq) {
	        nhpr = nlsq;
        }
        nren = 6;
        for(i=0; i<nren; i++) {
            ener[i]  =0.;
            enert[i] =0.;
            enert2[i]=0.;
        }

/* calculate and print initial temperature */

        eold = 0;
        for(j3=0; j3<nat3; j3++) {
	        eold += v[j3]*v[j3];
	    }
        eold *= (amas/2.);
        ener[1] = eold;
        if(ntt > 0) {
	        if(eold < 1.e-6) {
	            ekg = ekin0;
	        } else {
                ekg = eold;
	        }
        }

        double pres0=pres; // xxx ???

        *temp = ener[1]/fac;
        cout << "\n\n INITIAL TEMPERATURE IS :\n\n " << *temp << "\n";
/* dynamics step */
        for(nstep=0; nstep<nstlim; nstep++) {

/* put atoms in central periodic box */
            shia(nat, x, box, xmin);


            
/* calculate forces, potential energy, virial
 * and contribution to the radial distribution function
 */

            // 6.6.2 xxx
            // TODO: add asserts in main.cpp for input params

            // double gamma=1;

            if (ntt==3)
                langevin(nat, x, box, epslj, siglj, rcutf, &epot, f, &vir, rcutg, ngr, igr,
                gamma, amas, boltz, temp0, ig, v);
            else   
                forcea(nat, x, box, epslj, siglj, rcutf, &epot, f, &vir, rcutg, ngr, igr);


            ener[2] = epot;
            ener[3] = vir;

/* determine velocity scaling factor, when coupling to a bath */
            scal = (ntt == 1) ? sqrt(one+dtt*(ekin0/ekg-one)) : one;



            // 6.6.3 xxx 

            // double taup=1;
            // double betat=1;


            if (ntt==4)
                scal = pow((1 + betat*dt/taup * (pres - pres0)), 1/3);


/* perform leap-frog integration step,
 * calculate kinetic energy at time t-dt/2 and at time t,
 * and calculate pressure
 */
            eold = 0.;
            enew = 0.;
            for(j3=0; j3<nat3; j3++) {
                vh     = v[j3];
                vn     = (vh+f[j3]*dtm)*scal;
                eold  += vn*vn;
                enew  += (vh+vn)*(vh+vn);
                v[j3]  = vn;
                x[j3] += vn*dt;
            }
            eold *= (amas/2.);
            enew *= (amas/8.);
            ener[1] = enew;
            ener[0] = ener[1] + ener[2];
            pres = 2.*(enew-vir)/(vol*3.);
            ener[4] = pres;
            ener[5] = scal;
            if(ntt > 0){
		        ekg = eold;
		    }



            // 6.6.1 xxxx
            
            if (ntt==2){

                // time since last collision
                double tSinceColl = dt*nstep - lastcoll;
                // double sd = sqrt(boltz*temp0/amas)/3; // sd war zu gross, darum /3
                double sd = dtcoll/4;
                double t_rand = gauss(dtcoll, sd, &ig);


                // collide
                if (t_rand < tSinceColl){
                    // cout << "\n**********************\n";
                    // cout << "collision at step " << nstep << " \n";

                    int numOfColls = 1; // to implement (input)

                    for (int i=0; i<numOfColls; i++){

                        // choose atom
                        double r = a_random(&ig);
                        int atomN = (int) floor(r*nat);
                        int koord = atomN*3;

                        // randomize v[x/y/z]
                        double vGauss;
                        double mean = sqrt(3*boltz*temp0/amas);
                        sd = sqrt(boltz*temp0/amas);

                        vGauss = gauss(0, sd, &ig);
                        v[koord] = vGauss;

                        vGauss = gauss(0, sd, &ig);
                        v[koord+1] = vGauss;

                        vGauss = gauss(0, sd, &ig);
                        v[koord+2] = vGauss;
                    }

                    collCounter++;
                    lastcoll = dt*nstep;
                }
            }


		
/* update arrays for averages and fluctuations */
            t = t+dt;
            for(m=0; m<nren; m++) {
                enert[m]  += ener[m];
                enert2[m] += ener[m]*ener[m];
            }

/* write output per ... md steps */
		    if((nstep+1) == (nstep+1)/ntwx*ntwx && nstep < ntwxm){
		        packa (x, nat3, 1, trajectory_filename, ntpw);
		    }
		
            if(nstep == (nstep+1)/nhpr*nhpr)  {
		        cout.setf(ios::right);
		        cout << "\n\n";
                cout << setw(15) << "STEP";
                cout << setw(15) << "TIME";
                cout << setw(15) << "E-TOTAL";
                cout << setw(15) << "E-KINETIC";
                cout << setw(15) << "E-POTENTIAL";
                cout << setw(15) << "VIRIAL";
                cout << setw(15) << "PRESSURE";
                cout << setw(15) << "SCALE-T";
                cout << setw(15) << "TEMP";
                cout << "\n\n";
		    }
		    
            if((nstep+1) == (nstep+1)/ntpr*ntpr || nstep == 0) {
                cout << setw(15) << nstep+1 ;
                cout << setw(15) << t;
                for(k=0; k<nren; k++){
			        cout << setw(15) << ener[k] ;
			    }
                cout <<"\n";
            }


/* calculate and print center of mass motion
 * once in nlsq steps, at time t-dt/2
 */
            if((nstep+1) == (nstep+1)/nlsq*nlsq) {
		        for(j3=0; j3<nat3; j3++) {
		            f[j3] = x[j3]-v[j3]*dt5;
		        }
                cenma(nat, f, v, amas, &ekcm, xcm, vcm, 2);
            }
        }
/*
 * 27 mar 91
 * 4 or so lines down, there is a call to fabs(). Used to be abs().
 */
/* print averages */
        tspan = nstlim;
        for(m=0; m<nren; m++) {
            enert[m]  = enert[m]/tspan;
            enert2[m] = sqrt((double) fabs(enert2[m]/tspan-enert[m]*enert[m]));
        }
        cout << "\n\n AVERAGES ARE :\n\n";
        cout << "\n\n";
        cout.setf(ios::right);
        cout << setw(15) << "STEP";
        cout << setw(15) << "TIME";
        cout << setw(15) << "E-TOTAL";
        cout << setw(15) << "E-KINETIC";
        cout << setw(15) << "E-POTENTIAL";
        cout << setw(15) << "VIRIAL";
        cout << setw(15) << "PRESSURE";
        cout << setw(15) << "SCALE-T";
        cout << setw(15) << "Current Temp";
	    cout << "\n\n";	 
	
        cout.setf(ios::right);
        cout << setw(15) <<  nstlim;
        cout << setw(15) <<  t;
        for(k=0; k<nren; k++) {
	        cout << setw(15) <<  enert[k];
	    }

        cout << "\n";
        cout << "\n\n ROOT MEAN SQUARE FLUCTUATIONS :\n\n";
        cout << "\n\n";
        cout.setf(ios::right);
        cout << setw(15) << "STEP";
        cout << setw(15) << "TIME";
        cout << setw(15) << "E-TOTAL";
        cout << setw(15) << "E-KINETIC";
        cout << setw(15) << "E-POTENTIAL";
        cout << setw(15) << "VIRIAL";
        cout << setw(15) << "PRESSURE";
        cout << setw(15) << "SCALE-T";
        cout << "\n\n";
	
        cout.setf(ios::right);
        cout << setw(15) << nstlim;
        cout << setw(15) << t;
        for(k=0; k<nren; k++) {
	        cout << setw(15) << enert2[k];
	    }
        cout << setw(15) << ener[1];
	
        cout << "\n";
        *temp  = enert[1]/fac;
        temp2 = enert2[1]/fac;
        cout << "\n\n AVERAGE TEMPERATURE IS :\n\n " <<  *temp << "\n";
        cout << "\n\n ROOT MEAN SQUARE FLUCTUATIONS :\n\n " << temp2 << "\n";
        cout << "\n\n NUMBER OF COLLISIONS :\n\n " << collCounter << "\n"; // xxx
}


