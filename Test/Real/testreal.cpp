#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
//#include <complex>
#include <stdlib.h>
using namespace std;

#define LX 4
#define LY 4
#define D 2
#define PI 3.141592653589793
#define TWO_PI 6.283185307179586

#ifdef GPU
#include <accelmath.h> //For PGI
#endif

//typedef complex<double> Complex;
//#define cUnit Complex(1.0,0)


int main(int argc, char **argv) {

   double gaugeFree[LX][LY][D];

#pragma acc init

#pragma acc data copyout(gaugeFree[0:LX][0:LY][0:2])
 {
#pragma acc parallel loop collapse(2)
      for(int x=0; x<LX; x++)
          for(int y=0; y<LY; y++) {
              gaugeFree[x][y][0] = 1.0;
              gaugeFree[x][y][1] = 1.0;
               }
 } //acc pragma

  cout << "This code has run!" << endl;
  return 0;
} //main

