
#include <NTL/new.h>
#include "CC.h"

NTL_START_IMPL;

/* Implementation of complex arithmatic.
 *
 * Written by: Chris Studholme
 * Copyright:  GPL (http://www.fsf.org/copyleft/gpl.html)
 */


NTL_mod_polynomial_impl(RR,RRX,RRXModulus);

NTL_polymod_impl_noinit(RR,RRX,RRXModulus,CC);
NTL_eq_polymod_impl(RR,RRX,CC);
NTL_io_polymod_impl(RR,RRX,CC);
NTL_add_polymod_impl(RR,RRX,CC);
NTL_mul_polymod_impl(RR,RRX,CC);


RRXModulus* CC::F;
CC* CC::zero_ref;
class CC_initializer {
public:
  CC_initializer() {
    RRX f;
    SetCoeff(f,2);
    SetCoeff(f,0);
    CC::F = new RRXModulus(f);  // f = X^2+1
    CC::zero_ref = new CC();
  }
};
CC_initializer CC_initialize;


void sqr(CC& x, const CC& a) {
  RR u,v;
  sqr(u,real(a));
  sqr(v,imag(a));
  sub(u,u,v);
  mul(v,real(a),imag(a));
  add(v,v,v);
  x=u;
  SetImag(x,v);
}

void mul(CC& x, const CC& a, const CC& b) {
  // (ar+ai*i)*(br+bi*i) = (ar*br - ai*bi) + (ar*bi + ai*br)*i
  RR t,u,v;
  mul(u,real(a),real(b));
  mul(t,imag(a),imag(b));
  sub(u,u,t);
  mul(v,real(a),imag(b));
  mul(t,imag(a),real(b));
  add(v,v,t);
  x=u;
  SetImag(x,v);
}

CC operator/(const CC& a, const RR& b) {
  CC aa(a);
  aa.SetReal(real(a)/b);
  aa.SetImag(imag(a)/b);
  return aa;
}
CC operator/(const CC& a, const CC& b) {
  RR ar,ai;
  abs(ar,real(b));
  abs(ai,imag(b));
  RR nr,ni;
  RR t,d;
  if (ar<=ai) {
    t = real(b)/imag(b);
    d = imag(b)*(1+t*t);
    nr = (real(a)*t+imag(a))/d;
    ni = (imag(a)*t-real(a))/d;
  }
  else {
    t = imag(b)/real(b);
    d = real(b)*(1+t*t);
    nr = (real(a)+imag(a)*t)/d;
    ni = (imag(a)-real(a)*t)/d;
  }
  return CC(nr,ni);
}


void norm(RR& x, const CC& a) {
  RR t;
  sqr(x,real(a));
  sqr(t,imag(a));
  add(x,x,t);
}

void power(CC& x, const CC& a, long e) {
  if (e<=1) {
    if (e<0)
      Error("power: negative exponent");
    if (e==0)
      set(x);
    else if (&x!=&a) 
      x=a;
    return;
  }
  if (deg(rep(a))<=0) {
    x=power(real(a),e); // real only
    return;
  }

  CC aa(a);
  set(x);
  for (long i=NumBits(e)-1; i>=0; --i) { 
    sqr(x,x); 
    if (bit(e,i)) 
      mul(x,x,aa); 
  } 
}


/* This code from gcc file g++-3/std/complext.cc.
 */
void SqrRoot(CC& x, const CC& a) {
  if (IsZero(a)) {
    clear(x);
    return;
  }

  RR nr,ni;
  if (real(a)>0) {
    nr = sqrt(0.5*(abs(a)+real(a)));
    ni = imag(a)/nr/2;
  }
  else {
    ni = sqrt(0.5*(abs(a)-real(a)));
    if (imag(a)< 0)
      ni=-ni;
    nr=imag(a)/ni/2;
  }
  x=nr;
  x.SetImag(ni);
}



NTL_END_IMPL;
