#ifndef NTL_RRX__H
#define NTL_RRX__H

#include <NTL/RR.h>
//#include <NTL/polynomial.h>
#include "polynomial.h"
#include <NTL/vector.h>

NTL_OPEN_NNS;

extern const RR RR_zero;

NTL_vector_decl(RR,vec_RR);
NTL_eq_vector_decl(RR,vec_RR);
NTL_io_vector_decl(RR,vec_RR);

NTL_polynomial_decl(RR,vec_RR,RRX,RR_zero);
NTL_eq_polynomial_decl(RR,vec_RR,RRX);
NTL_io_polynomial_decl(RR,vec_RR,RRX);
NTL_add_polynomial_decl(RR,vec_RR,RRX);
NTL_mul_polynomial_decl(RR,vec_RR,RRX);
NTL_diff_polynomial_decl(RR,vec_RR,RRX);

inline void SetCoeff(RRX& x, long i, double d) {
  SetCoeff(x,i,to_RR(d));
}
inline void conv(RRX& x, double d) {
  conv(x,to_RR(d));
}

/* polynomial evaluation: promotion from double */
inline void eval(RR& x, const RRX& f, double d) {
  eval(x,f,to_RR(d));
}
inline RR eval(const RRX& f, double d) {
  RR x; eval(x,f,to_RR(d)); NTL_OPT_RETURN(RR,x); 
}


NTL_CLOSE_NNS;


#endif
