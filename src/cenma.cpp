/* w.f. van gunsteren, groningen, feb. 1988
 * h.j.w spoelder, Amsterdam, C-version 1.0, feb. 1988
 *
 *     cenma(nat, x, v, amas, ekcm, xcm, vcm, icm)
 *
 *     cenma will supply the center of mass coordinates (icm<>0).
 *     it also calculates the center of mass velocity and
 *     translational kinetic energy (icm<-1 or icm>+1).
 *
 *     nat = number of atoms
 *     x(1..3*nat) = atom cartesian coordinates
 *     v(1..3*nat) = atom velocities
 *     amas = mass of an atom
 *     ekcm = delivered with center of mass translational kinetic
 *            energy (if abs(icm)>1)
 *     xcm(1..3) = delivered with center of mass coordinates
 *                 (if icm<>0)
 *     vcm(1..3) = delivered with center of mass velocity
 *                 (if abs(icm)>1)
 *     icm =  0    : nothing is done
 *         = -1,+1 : xcm is calculated
 *         = -2,+2 : in addition, ekcm and vcm are calculated
 *                 : if icm>0, the results are printed
 *
 *
 *
 */
#include <iostream>      // For terminal I/O
#include <iomanip>       // For terminal output format
#include <cstdlib>
using namespace std;

void cenma( int nat, double x[], double v[], double amas, double *ekcm,
            double xcm[], double vcm[], int icm){

        double  tmas;
        int     j, m;

/*
 * initialise
 */
        if(!icm){
	  return;
	}
	
        tmas = nat*amas;
/*
 * calculate the center of mass coordinates
 */
        for(m=0; m<3; m++){
	  xcm[m] = 0.;
	}
	
        for(j=0; j<3*nat; j+=3) {
                for(m=0; m<3; m++)
                        xcm[m] += amas*x[j+m];
        }
        for(m=0; m<3; m++){
	  xcm[m] /= tmas;
	}
	
        if(icm > 0) {
                cout <<"\n\n   TOTAL MASS                       :    ";
                cout << setw(15) << setiosflags(ios::left) <<tmas ;
                cout << endl;
	        cout <<"   COORDINATES OF THE CENTER OF MASS:    ";
                for(m=0; m<3; m++)
		  cout << setw(15) << xcm[m];
                cout << endl;
        }
        if(abs(icm) < 2){
	  return;
        }

/*
 * calculate the velocity and the translational kinetic energy
 * of the center of mass
 */
        *ekcm = 0;
        for(m=0; m<3; m++){
	  vcm[m] = 0.;
	}
	
        for(j=0; j<3*nat; j+=3) {
                for(m=0; m<3; m++)
                        vcm[m] += amas*v[j+m];
        }
        for(m=0; m<3; m++) {
                vcm[m] /= tmas;
                *ekcm += vcm[m]*vcm[m];
        }
        *ekcm *= (tmas/2.);
        if(icm > 1) {
	        cout << "   VELOCITY OF THE CENTER OF MASS   :    ";
                for(m=0; m<3; m++){
		  cout << setw(15) << vcm[m];
		}
                cout << endl;
                cout << "   TRANSLATIONAL KINETIC ENERGY     :    " 
                     << setw(15) << setiosflags(ios::left) << *ekcm << "\n";
        }

      return;
}
