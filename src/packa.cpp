/*---------------------------FILE PACKA.C  #6 of 10-----------------------*/

/* w.f. van gunsteren, groningen, feb. 1988
 * h.j.w. spoelder, Amsterdam, C-version 1.0, feb. 1988
 * gee, C++ version 1.0, Oct. '01
 * mph, eth zurich, 01/2015, removed tape* as filename
 *
 * Note:
 *   Function coded below is substantially different from original 
 *   FORTRAN and C versions. It has been slimmed down to the minimum
 *   required by the rest of the MDATOM code. That is: 
 *        - binary write to file
 *        - formatted write to file
 *   (10/01, gee)
 *   (01/15, mph)
 *
 *     packa(r, ntot, mode, filename, ntpw)
 *
 *     packa writes or reads the array r to or from disk.
 *     it may write or read in binary form, in formatted form or may
 *     read in standard formatted form.
 *
 *     r(1.. ) = array to be written or read
 *     ntot = number of array elements in r to be written or read
 *            (ntpw<3)
 *          = not used (ntpw>2)
 *     mode = +1 : r will be written (ntpw<3)
 *          = -1 : r will be read (ntpw<4)
 *     filename = std::string
 *     ntpw = 0 : reading or writing in binary form
 *          = 1 : this option has been disabled
 *          = 2 : reading or writing in formatted form (10f8.4)
 *          = 3 : reading in standard formatted form    (%8lf %8lf %8lf)
 *
 *
 *
 */
#include  <fstream>   // For input/output file handling
#include  <iomanip>   // For I/O format
#include  <string>    // For filename storage
#include  <sstream>   // For appending int to string
                      // (filename formatting)
#include <cstdlib>                    
#include "mdprog.h"   // For run-time constants 
                      // and binary I/O function prototypes 

using namespace std;


void packa( double r[], int ntot, int mode, string filename, int ntpw){

  
	int     i, j;
        ofstream fileBW;          // Binary write
        ofstream fileFW;          // Formatted write
        
        if(ntpw == 1 || ntpw > 2) {
	     cerr << "Variable ntpw may only take values 0 or 2 ";
	     cerr << "in this program\n";
         return;
        }

        // Write out trajectory in binary form
        if(ntpw == 0) { 

	    fileBW.open(filename.c_str(),
                        ios::out | ios::app |
                        ios::binary);
            if(fileBW.bad()) {
	        cout << "I/O ERROR: cannot write to file: " 
                     << filename << endl;
                exit(1);
            }
            writedbl(fileBW, r, ntot);
            fileBW.close();
        }

    
        // Write out trajectory in formatted form
        if(ntpw == 2){

                fileFW.open(filename.c_str(), ios::out | ios::app);
                if(fileFW.bad()) {
	          cerr << "I/O ERROR: cannot write to file: " 
                       << filename << endl;
                  exit(1);
                }
                for(i=0; i<ntot; i += 10) {
		  for(j=i; (j<i+10 && j<ntot); j++){
		    fileFW << setw(10) << r[j];
		  }
                  fileFW << endl;
                }    
        }

}

