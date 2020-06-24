#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
using namespace std;

#define LX 512
#define LY 512
#define D 2
#define TWO_PI 6.283185307179586

#ifdef GPU
#include <accelmath.h>
#endif

#include "utils.h"

int main(int argc, char **argv) {

  double gauge[LX][LY][D][2];
  double psi1[LX][LY][D][2];
  double psi2[LX][LY][D][2];
  double mass = 0.0;
  double beta = 6.0;
  int numruns = 250000;

  zeroLat(gauge);
  zeroLat(psi1);
  gaussStart(gauge, beta);
  gaussComplex_F(psi1);
  double time0 = -((double)clock());

#pragma acc init

#pragma acc data copyout(psi1[0:LX][0:LY][0:D][0:2], gauge[0:LX][0:LY][0:D][0:2])
{

  for (int k=0; k<numruns; k++) {
    Dpsi(psi2, psi1, gauge, mass);
  }
 } //acc pragma
  double timeend = time0 + clock();
  cout << "Time for Mvps: " << timeend/CLOCKS_PER_SEC << endl;
 

}

  




