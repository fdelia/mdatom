#include <iostream>
#include "mdprog.h"
using namespace std;

//
// Auxiliary functions for program <mdatom>.
// 
// Function list:
//     readdbl  - called from main.cpp and packa.cpp
//     writedbl - called from main.cpp and packa.cpp
//     usage    - called from main.cpp
//

//-------------------------------------------------
// readdbl - read binary input stream as data of type
//           double into an array
//         - such a function not present in std library
// 
// Parameter list:
//     in    - input stream
//     value - target array
//     length- nr. of cells in target array
//-------------------------------------------------
void readdbl(istream& in, double value[], int length)
{
  in.read(reinterpret_cast<char*>(value), length*sizeof(double));
}

//-------------------------------------------------
// writedbl - write contents of an array of data of type
//           double to output stream as binary data.
//          - such a function not present in standard library
// 
// Parameter list:
//     out    - output stream
//     value - target array
//     length- nr. of cells in target array
//-------------------------------------------------
void writedbl(ostream& out, double value[], int length)
{
  out.write(reinterpret_cast<char*>(value), length*sizeof(double));
}

//-------------------------------------------------
// usage - write out usage message
//-------------------------------------------------
void usage ()
{
      cerr << "Usage promda input coord > output \n";
}


//-------------------------------------------------
// space - write out white space
//-------------------------------------------------
void space(int n) {
        while(n--)
                cout << " ";
}
        
