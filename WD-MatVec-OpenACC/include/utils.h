#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

// Zero lattice 2D
template<typename T> inline void zeroLat(T v[LX][LY][2][2]) {
  for(int x=0; x<LX; x++) {
    for(int y=0; y<LY; y++) {
      for(int mu=0; mu<2; mu++) {
  v[x][y][mu][0] = 0.0;
  v[x][y][mu][1] = 0.0;
  }
 }
}
}

void gaussComplex_F(double field[LX][LY][2][2]) {

  double r, theta, sum;
  for(int x=0; x<LX; x++) {
    for(int y=0; y<LY; y++){
      r = sqrt(-2.0*log((double(rand())/ (RAND_MAX))));
      theta = TWO_PI*(double(rand())/ (RAND_MAX));
      field[x][y][0][0] = r*cos(theta);
      field[x][y][0][1] = r*sin(theta);

      r = sqrt(-2.0*log((double(rand())/ (RAND_MAX))));
      theta = TWO_PI*(double(rand())/ (RAND_MAX));
      field[x][y][1][0] = r*cos(theta);
      field[x][y][1][1] = r*sin(theta);
    }
 }
    return;
}

void gaussStart(double gauge[LX][LY][2][2], double beta){

   for(int x=0; x<LX; x++) {
     for(int y=0; y<LY; y++) {
	for (int mu=0; mu<2; mu++) {
       gauge[x][y][mu][0] = cos(sqrt(1.0/beta))*(double(rand())/ (RAND_MAX));
       gauge[x][y][mu][1] = sin(sqrt(1.0/beta))*(double(rand())/ (RAND_MAX));
     }
   }
  }
   return;
}

// Wilson stencil
// // D_{W}(n,m) = (m_{0} + 2r)\delta(n,m)
// //               - 1/2 Sum [(r-\sigma_{mu}) U_{n,\mu} \delta_{n,m-\hat{\mu}} +
// //                          (r+\sigma_{mu}) U^{\dagger}_{m,\mu} \delta_{n,m+\hat{\mu}}]
// //
// // sigma_1 = | 0  1 |  sigma_2 = | 0 -i | sigma_3 = i*sigma_1*sigma_2 = | 1  0 |
// //           | 1  0 |            | i  0 |                               | 0 -1 |
//
//

void Dpsi(double psi2[LX][LY][2][2], double psi1[LX][LY][2][2],
     double gauge[LX][LY][2][2], double mass){

  double m0 = mass;
  double  r = 1.0;
  double constant = (2*r + m0);
  int xp1, xm1, yp1, ym1;
#pragma acc data present(psi1[0:LX][0:LY][0:D][0:2], gauge[0:LX][0:LY][0:D][0:2]) create(psi2[0:LX][0:LY][0:D][0:2])
 {

#pragma acc parallel loop tile(32,2)
  for(int x=0; x<LX; x++) { 
    xp1 = (x+1)%LX;
    xm1 = (x-1+LX)%LX;
    for(int y=0; y<LY; y++) {
      yp1 = (y+1)%LY;
      ym1 = (y-1+LY)%LY;

	//upper
	psi2[x][y][0][0] = constant * psi1[x][y][0][0] -
	
	0.5*( r * (gauge[x][y][0][0] * psi1[x][y][0][0] - gauge[x][y][0][1] * psi1[x][y][0][1]) - 
	    
	    (gauge[x][y][0][0] * psi1[xp1][y][1][0] - gauge[x][y][0][1] * psi1[xp1][y][1][1]) +

	    ( r * (gauge[xp1][y][0][0] * psi1[x][y][0][0] + gauge[xp1][y][0][1] * psi1[x][y][0][1])) +

            (gauge[xp1][y][0][0] * psi1[xp1][y][1][0] + gauge[xp1][y][0][1] * psi1[xp1][y][1][1]) +

	    ( r * (gauge[x][y][1][0] * psi1[x][yp1][0][0] - gauge[x][y][1][1] * psi1[x][yp1][0][1])) +

           ( r * (gauge[x][yp1][1][0] * psi1[x][yp1][0][0] + gauge[x][yp1][1][1] * psi1[x][yp1][0][1])));

	psi2[x][y][0][1] = constant * psi1[x][y][0][1] -
	
	0.5*( r * (gauge[x][y][0][0] * psi1[x][y][0][1] + gauge[x][y][0][1] * psi1[x][y][0][0]) -

	    (gauge[x][y][0][0] * psi1[xp1][y][1][1] + gauge[x][y][0][1] * psi1[xp1][y][1][0]) +

	    ( r * (gauge[xp1][y][0][0] * psi1[x][y][0][0] - gauge[xp1][y][0][1] * psi1[x][y][0][1])) +

            (gauge[xp1][y][0][0] * psi1[xp1][y][1][0] - gauge[xp1][y][0][1] * psi1[xp1][y][1][1]) +

	    ( r * (gauge[x][y][1][0] * psi1[x][yp1][0][1] + gauge[x][y][1][1] * psi1[x][yp1][0][0])) +

	    ( r * (gauge[x][yp1][1][0] * psi1[x][yp1][0][1] - gauge[x][yp1][1][1] * psi1[x][yp1][0][0])) +

	    (gauge[x][y][1][0] * psi1[x][yp1][1][0] + gauge[x][y][1][1] * psi1[x][yp1][1][1]) +

	    (gauge[x][y][1][0] * psi1[x][yp1][1][1] - gauge[x][y][1][1] * psi1[x][yp1][1][0]));

	//Lower
	psi2[x][y][1][0] = constant * psi1[x][y][1][0] -

	0.5*( r * (gauge[x][y][0][0] * psi1[xp1][y][1][0] - gauge[x][y][0][1] * psi1[xp1][y][1][1]) -

	    (gauge[x][y][0][0] * psi1[xp1][y][0][0] - gauge[x][y][0][1] * psi1[xp1][y][0][1]) +
	 
	    (gauge[xm1][y][0][0] * psi1[xm1][y][0][0] + gauge[xm1][y][0][1] * psi1[xm1][y][0][1]) +

	    (r * (gauge[xm1][y][0][0] * psi1[xm1][y][1][0] + gauge[xm1][y][0][1] * psi1[xm1][y][1][1])) +

	    ( r * (gauge[x][y][1][0] * psi1[x][yp1][1][0] - gauge[x][y][1][1] * psi1[x][yp1][1][1])) +

	    ( r * (gauge[x][ym1][1][0] * psi1[x][ym1][1][0] + gauge[x][ym1][1][1] * psi1[x][ym1][1][1])));

	psi2[x][y][1][1] = constant * psi1[x][y][1][1] - 

	0.5*( r * (gauge[x][y][0][0] * psi1[xp1][y][1][1] + gauge[x][y][0][1] * psi1[xp1][y][1][0]) -

	    (gauge[x][y][0][0] * psi1[xp1][y][0][1] + gauge[x][y][0][1] * psi1[xp1][y][0][1]) +

	    (gauge[xm1][y][0][0] * psi1[xm1][y][0][1] - gauge[xm1][y][0][1] * psi1[xm1][y][0][0]) +

	    (r * (gauge[xm1][y][0][0] * psi1[xm1][y][1][1] - gauge[xm1][y][0][1] * psi1[xm1][y][1][0])) +

	    ( r * (gauge[x][y][1][0] * psi1[x][yp1][1][1] + gauge[x][y][1][1] * psi1[x][yp1][1][0])) +

	    ( r * (gauge[x][ym1][1][0] * psi1[x][ym1][1][1] - gauge[x][ym1][1][1] * psi1[x][ym1][1][0])) -

	    (gauge[x][y][1][0] * psi1[x][yp1][0][0] + gauge[x][y][1][1] * psi1[x][yp1][1][0]) +

	    (gauge[x][ym1][1][0] * psi1[x][ym1][0][1] - gauge[x][ym1][1][1] * psi1[x][ym1][0][1]));

	    
   } //y
 } //x
} //acc pragma
} //void

#endif
