#ifndef NTL_CC__H
#define NTL_CC__H

#include "RRX.h"
#include "polymod.h"

/*
Copyright (C) 2004 Chris Studholme

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/


NTL_OPEN_NNS;

NTL_mod_polynomial_decl(RR,RRX,RRXModulus);

NTL_polymod_decl_begin(RR,RRX,RRXModulus,CC);
  CC(const RR& real, const RR& imag) {
    SetImag(imag);
    SetReal(real);
  }

  // set parts
  inline void SetReal(const RR& a) {
    SetCoeff(rep,0,a); 
  }
  inline void SetImag(const RR& a) {
    SetCoeff(rep,1,a); 
  }

  // get parts (read-only reference)
  inline const RR& CC::real() {
    return rep[0];
  }
  inline const RR& CC::imag() {
    return rep[1];
  }

  // read-only reference to i
  static const CC& CC::I() { 
    static CC i(to_RR(0),to_RR(1));    
    return i;
  }
NTL_polymod_decl_end(RR,RRX,RRXModulus,CC);

NTL_eq_polymod_decl(RR,RRX,CC);
NTL_io_polymod_decl(RR,RRX,CC);
NTL_add_polymod_decl(RR,RRX,CC);
NTL_mul_polymod_decl(RR,RRX,CC);

// set parts
inline void SetReal(CC& x, const RR& a) {
  x.SetReal(a);
}
inline void SetImag(CC& x, const RR& a) {
  x.SetImag(a);
}

// real part of x
inline const RR& real(const CC& x) {
  return rep(x)[0];
}

// imaginary part of x
inline const RR& imag(const CC& x) {
  return rep(x)[1];
}

// squaring
void sqr(CC& x, const CC& a);
inline CC sqr(const CC& a) {
  CC x;  sqr(x,a);  NTL_OPT_RETURN(CC,x);
}

// multiplication
void mul(CC& x, const CC& a, const CC& b);
inline CC operator*(const CC& a, const CC& b) {
  CC x;  mul(x,a,b);  NTL_OPT_RETURN(CC,x);
}
inline CC operator*=(CC& x, const CC& a) {
  mul(x,x,a);
  return x;
}


// division
CC operator/(const CC& a, const RR& b);
CC operator/(const CC& a, const CC& b);



// square root
void SqrRoot(CC& x, const CC& a);
inline CC SqrRoot(const CC& a) {
  CC x;  SqrRoot(x,a);  NTL_OPT_RETURN(CC,x);
}
inline CC sqrt(const CC& a) {
  CC x;  SqrRoot(x,a);  NTL_OPT_RETURN(CC,x);
}


// norm and absolute value
void norm(RR& x, const CC& a);
inline RR norm(const CC& a) {
  RR x;  norm(x,a);  NTL_OPT_RETURN(RR,x);
}
inline void abs(RR& x, const CC& a) {
  SqrRoot(x,norm(a));
}
inline RR abs(const CC& a) {
  RR x;  abs(x,a);  NTL_OPT_RETURN(RR,x);
}


inline void polar(CC& x, const RR& r, const RR& t) {
  RR i;
  mul(i,r,sin(t));
  x = r*cos(t);
  x.SetImag(i);
}
inline CC polar(RR r, RR t) {
  CC x; polar(x,r,t); NTL_OPT_RETURN(CC,x); 
}

// exp
inline void exp(CC& x, const CC& a) {
  polar(x,exp(real(a)),imag(a));
}
inline CC exp(const CC& a) {
  CC x;  exp(x,a);  NTL_OPT_RETURN(CC,x);
}


// power
void power(CC& x, const CC& a, long e);
inline CC power(const CC& a, long e) {
  CC x;  power(x,a,e);  NTL_OPT_RETURN(CC,x);
}



NTL_CLOSE_NNS;

#endif
