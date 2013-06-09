#ifndef NTL_polymod__H
#define NTL_polymod__H

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

*****************************************************

MODULE: polymod

SUMMARY:

Macros are defined providing template-like classes for polynomials modulo
a fixed polynomial modulus.

  ... see polymod.txt for complete documentation ...

*/


#include <NTL/tools.h>


#define NTL_mod_polynomial_decl(T,TX,TXModulus) \
class TXModulus { \
public: \
  /* initially in an unusable state */ \
  TXModulus() {}  \
 \
  /* copy constructor */ \
  TXModulus(const TXModulus& other) { \
    f=other.f; \
  } \
 \
  /* initialize with f, deg(f) > 0 */ \
  TXModulus(const TX& f); \
 \
  ~TXModulus() {} \
 \
  /* assignment operator */ \
  inline TXModulus& operator=(const TXModulus& other) { \
    f=other.f; \
    return *this; \
  } \
 \
  /* read-only access to f, implicit conversion operator */ \
  inline operator const TX&() const {  \
    return f; \
  } \
 \
  /* read-only access to f, explicit notation */ \
  inline const TX& val() const {  \
    return f; \
  } \
 \
public: \
  TX f; /* modulus */ \
 \
}; \
 \
/* pre-computes information about f and stores it in F. \
 * Note that the declaration TXModulus F(f) is equivalent to \
 * TXModulus F; build(F, f). \
 */ \
void build(TXModulus& F, const TX& f); \
 \
/* In the following, f refers to the polynomial f supplied to the \
 * build routine, and n = deg(f). \
 */ \
 \
/* return n=deg(f) */ \
inline long deg(const TXModulus& F) { \
  return deg(F.f); \
} \



#define NTL_mod_polynomial_impl(T,TX,TXModulus) \
TXModulus::TXModulus(const TX& f) { \
  if (deg(f)<=0) \
    Error("TXModulus: deg(f) must be at least 1"); \
  this->f = f; \
} \
void build(TXModulus& F, const TX& f) { \
  if (deg(f)<=0) \
    Error("build: deg(f) must be at least 1"); \
  F.f = f; \
} \



/* Declaration of polymod class.
 */
#define NTL_polymod_decl(T,TX,TXModulus,TE) \
NTL_polymod_decl_begin(T,TX,TXModulus,TE) \
NTL_polymod_decl_end(T,TX,TXModulus,TE) \



#define NTL_polymod_decl_begin(T,TX,TXModulus,TE) \
class TE { \
/*protected:*/ \
public: \
  static TXModulus* F; /* modulus */ \
  static TE* zero_ref; \
 \
  TX rep; /* polynomial residue */ \
 \
public: \
  /* Default constructor.  Initial value is 0. \
   */ \
  TE() {} \
 \
  /* Copy constructor.  \
   */ \
  TE(const TE& other) {  \
    rep=other.rep; \
  } \
 \
  /* Initialize to 0 and swap with other.  \
   */ \
  TE(TE& other, NTL_NNS INIT_TRANS_TYPE) { \
    NTL_NNS swap(rep,other.rep); \
  } \
 \
  /* Destructor.  Free all allocated memory. \
   */  \
  ~TE() {} \
 \
  /* Assignment from existing polynomial.   \
   */ \
  inline TE& operator=(const TE& other) { \
    rep=other.rep; \
    return *this; \
  } \
 \
  /* Assignment from constant term.   \
   */ \
  inline TE& operator=(const T& c) { \
    conv(rep,c); \
    return *this; \
  } \
 \
  /* Release all allocated space and set to 0. \
   */ \
  inline void kill() {  \
    rep.kill(); \
  } \
 \
  /* Set to zero. \
   */ \
  inline void clear() { \
    NTL_NNS clear(rep); \
  } \
 \
  /* Swaps this and other by swapping internal pointers.  \
   */ \
  inline void swap(TE& other) {  \
    NTL_NNS swap(rep,other.rep); \
  } \
 \
  /* Read-only reference to zero polynomial.  \
   */ \
  static inline const TE& zero() { \
    return *zero_ref; \
  } \
 \
  static void init(const TX& f); \
  static void init(const TXModulus& F); \
 \
  static inline const TXModulus& modulus() { \
    return *F; \
  } \
 \
  static inline long degree() { \
    return deg(*F); \
  } \
 \
public: \



#define NTL_polymod_decl_end(T,TX,TXModulus,TE) \
};\
 \
/* miscellany */ \
inline const TX& rep(const TE& a) { \
  return a.rep; \
} \
 \
inline bool IsZero(const TE &a) { \
  return IsZero(rep(a)); \
} \
 \
inline void clear(TE &x) { \
  x.clear(); \
} \
 \
inline void swap(TE &x, TE &y) { \
  x.swap(y); \
} \
 \
void conv(TE& dest, const TX& src); \
 \
inline void conv(TE& dest, const T& src) { \
  dest=src; \
} \



#define NTL_bak_polymod_decl(T,TX,TXModulus,TE,TEBak) \
class TEBak { \
public: \
  TEBak(); \
  ~TEBak(); \
  void save(); \
  void restore(); \
 \
private: \
  TXModulus* backup; \
  bool auto_restore; \
 \
  TEBak(const TEBak&) {}  /* disabled */ \
  void operator=(const TEBak&) {}  /* disabled */ \
}; \




#define NTL_eq_polymod_decl(T,TX,TE)  \
inline bool operator==(const TE& a, const TE& b) { \
  return a.rep==b.rep; \
} \
inline bool operator!=(const TE& a, const TE& b) { \
  return a.rep!=b.rep; \
} \
 \
/* promotion from T */ \
inline bool operator==(const T& a, const TE& b) { \
  return a==rep(b); \
} \
inline bool operator==(const TE& a, const T& b) { \
  return rep(a)==b; \
} \
inline bool operator!=(const T& a, const TE& b) { \
  return a!=rep(b); \
} \
inline bool operator!=(const TE& a, const T& b) { \
  return rep(a)!=b; \
} \



#define NTL_io_polymod_decl(T,TX,TE)  \
NTL_SNS istream& operator>>(NTL_SNS istream&, TE&); \
NTL_SNS ostream& operator<<(NTL_SNS ostream&, const TE&); \




// addition, subtraction and negation methods
#define NTL_add_polymod_decl(T,TX,TE)  \
inline void add(TE& x, const TE& a, const TE& b) { \
  add(x.rep,rep(a),rep(b)); \
} \
inline void sub(TE& x, const TE& a, const TE& b) { \
  sub(x.rep,rep(a),rep(b)); \
} \
inline void negate(TE& x, const TE& a) { \
  negate(x.rep,rep(a)); \
} \
 \
inline TE operator+(const TE& a, const TE& b) { \
  TE x; add(x,a,b); NTL_OPT_RETURN(TE, x); \
} \
inline TE operator-(const TE& a, const TE& b) { \
  TE x; sub(x,a,b); NTL_OPT_RETURN(TE, x); \
} \
inline TE operator-(const TE& a) { \
  TE x; negate(x,a); NTL_OPT_RETURN(TE, x); \
} \
 \
inline TE& operator+=(TE& x, const TE& a) { \
  add(x,x,a);  return x; \
} \
 \
inline TE& operator-=(TE& x, const TE& a) { \
  sub(x,x,a);  return x; \
} \
 \
/* promotions from T */ \
inline void add(TE& x, const T& a, const TE& b) { \
  add(x.rep,a,rep(b)); \
} \
inline void add(TE& x, const TE& a, const T& b) { \
  add(x.rep,rep(a),b); \
} \
inline void sub(TE& x, const T& a, const TE& b) { \
  sub(x.rep,a,rep(b)); \
} \
inline void sub(TE& x, const TE& a, const T& b) { \
  sub(x.rep,rep(a),b); \
} \
 \
inline TE operator+(const T& a, const TE& b) { \
  TE x; add(x,a,b); NTL_OPT_RETURN(TE,x); \
} \
inline TE operator+(const TE& a, const T& b) { \
  TE x; add(x,a,b); NTL_OPT_RETURN(TE,x); \
} \
inline TE operator-(const T& a, const TE& b) { \
  TE x; sub(x,a,b); NTL_OPT_RETURN(TE,x); \
} \
inline TE operator-(const TE& a, const T& b) { \
  TE x; sub(x,a,b); NTL_OPT_RETURN(TE,x); \
} \
 \
inline TE& operator+=(TE& x, const T& a) { \
  add(x,x,a); return x; \
} \
 \
inline TE& operator-=(TE& x, const T& a) { \
  sub(x,x,a); return x; \
} \



// scalar multiplication
#define NTL_mul_polymod_decl(T,TX,TE)  \
inline void set(TE& x) { \
  set(x.rep); \
} \
inline bool IsOne(const TE& x) { \
  return IsOne(rep(x)); \
} \
inline void SetX(TE& x) { \
  SetX(x.rep); \
} \
inline bool IsX(const TE& x) { \
  return IsX(rep(x)); \
} \
 \
inline void mul(TE& x, const TE& a, const T& b) { \
  mul(x.rep,rep(a),b); \
} \
inline void mul(TE& x, const T& a, const TE& b) { \
  mul(x.rep,a,rep(b)); \
} \
 \
inline TE operator*(const T& a, const TE& b) { \
  TE x; mul(x,a,b); NTL_OPT_RETURN(TE, x); \
} \
inline TE operator*(const TE& a, const T& b) { \
  TE x; mul(x,a,b); NTL_OPT_RETURN(TE, x); \
} \
inline TE operator*=(TE& x, const T& a) { \
  mul(x,x,a); \
  return x; \
} \




/**************** implementation macros ****************/


#define NTL_polymod_impl(T,TX,TXModulus,TE) \
NTL_polymod_impl_noinit(T,TX,TXModulus,TE) \
TXModulus* TE::F=NULL; \
TE* TE::zero_ref=NULL; \



#define NTL_polymod_impl_noinit(T,TX,TXModulus,TE) \
void TE::init(const TXModulus& F) { \
  if (TE::F) delete TE::F; \
  if (zero_ref) delete zero_ref; \
  TE::F = new TXModulus(F); \
  zero_ref = new TE(); \
} \
void TE::init(const TX& f) { \
  if (F) delete F; \
  if (zero_ref) delete zero_ref; \
  F = new TXModulus(f); \
  zero_ref = new TE(); \
} \
void conv(TE& x, const TX& a) { \
  const TX& f = TE::F->val(); \
  long n = deg(f); \
  if (deg(a)<n) { \
    x.rep = a; \
    return; \
  } \
  /* reduce a modulo f */ \
  trunc(x.rep,a,n); \
  /* fn is X^n mod f */ \
  TX fn; \
  trunc(fn,f,n); \
  negate(fn,fn); \
  /* add in a[n]*X^n mod f */ \
  TX tmp; \
  if (!IsZero(a[n])) { \
    mul(tmp,a[n],fn); \
    x.rep+=tmp; \
  } \
  TX fi(fn); \
  for (long i=n+1; i<=deg(a); ++i) { \
    /* fi is X^i mod f */ \
    fi<<=1; \
    mul(tmp,fi[n],fn); \
    fi+=tmp; \
    /* add in a[i]*X^i mod f */ \
    if (!IsZero(a[i])) { \
      mul(tmp,a[i],fi); \
      x.rep+=tmp; \
    } \
  } \
} \


#define NTL_bak_polymod_impl(T,TX,TXModulus,TE,TEBak) \
TEBak::TEBak() { \
  backup=NULL; \
  auto_restore=false; \
} \
TEBak::~TEBak() { \
  if (auto_restore) \
    restore(); \
  if (backup) delete backup; \
} \
void TEBak::save() { \
  if (backup) delete backup; \
  backup = TE::F ? new TXModulus(*TE::F) : NULL; \
  auto_restore = true; \
} \
void TEBak::restore() { \
  if (backup) \
    TE::init(*backup); \
  else { \
    if (TE::zero_ref) { \
      delete TE::zero_ref; \
      TE::zero_ref = NULL; \
    } \
    if (TE::F) { \
      delete TE::F; \
      TE::F = NULL; \
    } \
  } \
  auto_restore=false; \
} \



#define NTL_eq_polymod_impl(T,TX,TE) \
/* all methods inline, reserved for future use */ \



#define NTL_io_polymod_impl(T,TX,TE) \
NTL_SNS istream & operator>>(NTL_SNS istream& in, TE& a) { \
  TX x; \
  in>>x; \
  conv(a,x); \
  return in; \
} \
  \
NTL_SNS ostream& operator<<(NTL_SNS ostream& out, const TE& a) { \
  out<<a.rep; \
  return out; \
} \



#define NTL_add_polymod_impl(T,TX,TE)  \
/* all methods inline, reserved for future use */ \



// scalar multiplication
#define NTL_mul_polymod_impl_plain(T,TX,TE)  \
/* all methods inline, reserved for future use */ \



#define NTL_mul_polymod_impl(T,vec_T,TX) \
NTL_mul_polymod_impl_plain(T,vec_T,TX); \




#endif
