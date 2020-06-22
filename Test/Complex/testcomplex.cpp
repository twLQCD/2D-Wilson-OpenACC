#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>
#include <stdlib.h>
using namespace std;

#define LX 64
#define LY 64
#define D 2
#define PI 3.141592653589793
#define TWO_PI 6.283185307179586

#ifdef GPU
#include <accelmath.h> //For PGI
#endif

typedef complex<double> Complex;
#define cUnit Complex(1.0,0)


int main(int argc, char **argv) {

   Complex gaugeFree[LX][LY][D];

#pragma acc init

#pragma acc data copyout(gaugeFree[0:LX][0:LY][0:2])
 {
#pragma acc parallel loop collapse(2)
      for(int x=0; x<LX; x++)
          for(int y=0; y<LY; y++) { 
              gaugeFree[x][y][0] = cUnit;
               gaugeFree[x][y][1] = cUnit;
               }
 } //acc pragma
  return 0;
} //main  
