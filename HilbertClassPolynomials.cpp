#include <math.h>

//#include <std/complext.cc>

#include "RRX.h"
#include "CC.h"

#include "HilbertClassPolynomials.h"


/* Calculation of Hilbert class polynomials.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */


/**************** class HilbertClassPolynomials ****************/

const ZZX& HilbertClassPolynomials::poly(long D) {
  if (D>-3) {
    if (D<3) 
      Error("HilbertClassPolynomials: invalid discrimant");
    D=-D;
  }
  if (D&2)
    Error("HilbertClassPolynomials: invalid discrimant");
  ZZX& f = cached_poly(D);
  if (IsZero(f))
    get_poly(f,D);
  return f;
}

ZZX& HilbertClassPolynomials::cached_poly(long D) {
  long i=(D&1)+(-D-3)/2;
  if (i>=cache.length())
    cache.SetLength(i+1);
  return cache.RawGet(i);
}



/**************** class HCP_generate ****************/


inline CC function_D(const CC& tau) {
  RR pi2;
  ComputePi(pi2);
  pi2*=2;
  CC q;
  q = exp(CC::I()*tau*pi2);
  CC r,t;
  set(r);
  for (long n=1; n<100; ++n) {
    power(t,q,n*(3*n-1)/2);
    t += power(q,n*(3*n+1)/2);
    r += n&1 ? -t : t;
  }
  return q*power(r,24);
}

CC function_j(const CC& tau) {
  CC ftau;
  ftau = function_D(to_RR(2)*tau)/function_D(tau);
  CC result;
  result = power(to_RR(256)*ftau+to_RR(1),3)/ftau;
  return result;
}


/* HilbertClassPolynomial.  (Algorithm 7.6.1 in Cohen)
 * Requires: D negative discriminant.
 * Computes a monic polynomial P of degree h(D) of which j((D+D^{1/2})/2) is
 * a root.
 */
void HCP_generate::get_poly(ZZX& P, long D) {
  // step 1
  RRX Pr;
  set(Pr);
  long b = D&1;
  long B;
  B=(long)sqrt(-D/3.0);
  CC Dc;
  Dc = to_RR(D);

  ZZ t;
  RRX T;
  do {
    // step 2
    t = (b*b-D)/4;
    long a = max(b,1);
  
    do {
      // step 3
      if (divide(t,a)) {
	CC j;
	j = function_j((to_RR(-b)+sqrt(Dc))/to_RR(2*a));
	cout<<"HCP: a="<<a<<"  b="<<b<<"  j="<<j<<endl;
	if (a==b || b==0 || a*a==t) {
	  clear(T);
	  SetCoeff(T,0,-real(j));
	  SetCoeff(T,1);
	  //P = P * (X-j);  
	  Pr *= T;
	}
	else {
	  clear(T);
	  SetCoeff(T,0,norm(j));
	  SetCoeff(T,1,-2*real(j));
	  SetCoeff(T,2);
	  //P = P * (X^2 - 2*Re(j)*X + Norm(j));
	  Pr *= T;
	}
      }

      // step 4
      ++a;
    } while (a*a<=t);

    // step 5
    b+=2;
  } while (b<=B);
  
  // round coefficients of Pr
  P.SetMaxLength(deg(Pr)+1);
  clear(P);
  for (long i=0; i<=deg(Pr); ++i)
    SetCoeff(P,i,RoundToZZ(coeff(Pr,i)));
}
