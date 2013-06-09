#ifndef NTL_polynomial__H
#define NTL_polynomial__H

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

MODULE: polynomial

SUMMARY:

Macros are defined providing template-like classes for single variable
polynomials.

  ... see polynomial.txt for complete documentation ...

*/


#include <NTL/tools.h>



#ifndef NTL_RANGE_CHECK
#define NTL_POLY_RANGE_CHECK_CODE 
#else
#define NTL_POLY_RANGE_CHECK_CODE  if (i<0) RangeError(i)
#endif


/* Declaration of polynomial class.
 */
#define NTL_polynomial_decl(T,vec_T,TX,zero_T) \
NTL_polynomial_decl_begin(T,vec_T,TX,zero_T) \
NTL_polynomial_decl_end(T,vec_T,TX,zero_T) \



#define NTL_polynomial_decl_begin(T,vec_T,TX,zero_T) \
class TX { \
public: \
  vec_T rep; /* coefficients (should be protected too) */ \
 \
  /* Default constructor.  Initial value is 0. \
   */ \
  TX() {} \
 \
  /* Construct from single constant term.  Initial value is X^i*c.  \
   */ \
  TX(long i, const T& c) { \
    if (!IsZero(c)) { \
      rep.SetLength(i+1); \
      rep.RawGet(i)=c; \
      while (i>0) NTL_NNS clear(rep.RawGet(--i)); \
    } \
  } \
 \
  /* Copy constructor.  If normalize==true, leading zeros are stripped  \
   * from rep. \
   */ \
  TX(const TX& other, bool normalize=false);  \
 \
  /* Initialize to 0 and swap with other.  \
   */ \
  TX(TX& other, NTL_NNS INIT_TRANS_TYPE) { \
    NTL_NNS swap(rep,other.rep); \
  } \
 \
  /* Destructor.  Free all allocated memory. \
   */  \
  ~TX() {} \
 \
  /* Assignment from existing polynomial.   \
   */ \
  inline TX& operator=(const TX& other) { \
    rep=other.rep; \
    return *this; \
  } \
 \
  /* Assignment from constant term.   \
   */ \
  inline TX& operator=(const T& c) { \
    if (IsZero(c)) \
      rep.SetLength(0); \
    else { \
      /* c could be element of rep, but that's ok */ \
      rep.SetLength(1); \
      rep.RawGet(0)=c; \
    } \
    return *this; \
  } \
 \
  /* Strip leading zeros from rep.  \
   */  \
  inline void normalize() { \
    long i = rep.length(); \
    while (i>0 && IsZero(rep.RawGet(--i))) \
      rep.SetLength(i); \
  } \
 \
  /* Direct access to coefficients.  \
   */ \
  inline const T& operator[](long i) const { \
    return (0<=i && i<rep.length()) ? rep.RawGet(i) : zero_T; \
  } \
 \
  /* Release all allocated space and set to 0. \
   */ \
  inline void kill() {  \
    rep.kill(); \
  } \
 \
  /* Increase number of coordinates to n and zero new leading coordinates. \
   * Does nothing if n is less than current length. \
   */ \
  inline void ExtendTo(long n) { \
    long cur = rep.length(); \
    if (n>cur) { \
      rep.SetLength(n); \
      while (cur<n) NTL_NNS clear(rep.RawGet(cur++)); \
    } \
  } \
 \
  /* Pre-allocate spaces for n coefficients.  The polynomial represented  \
   * is not changed.  \
   */ \
  inline void SetMaxLength(long n) {  \
    rep.SetMaxLength(n); \
  } \
 \
  /* Maximum length ever achieved.  \
   */ \
  inline long MaxLength() const { \
    return rep.MaxLength(); \
  } \
 \
  /* Set to zero. \
   */ \
  inline void clear() { \
    rep.SetLength(0); \
  } \
 \
  /* Swaps this and other by swapping internal pointers.  \
   */ \
  inline void swap(TX& other) {  \
    NTL_NNS swap(rep,other.rep); \
  } \
 \
  /* Read-only reference to zero polynomial.  \
   */ \
  static inline const TX& zero() { \
    static TX zero_ref; \
    return zero_ref; \
  } \
public: \



#define NTL_polynomial_decl_end(T,vec_T,TX,zero_T) \
};\
 \
/* miscellany */ \
inline bool IsZero(const TX &x) { \
  return x.rep.length()==0; \
} \
 \
inline void clear(TX &x) { \
  x.clear(); \
} \
 \
inline void swap(TX &x, TX &y) { \
  swap(x.rep,y.rep); \
} \
 \
inline void conv(TX& dest, const T& src) { \
  dest=src; \
} \
 \
inline void conv(TX& dest, const vec_T& src) { \
  dest.rep=src; \
  dest.normalize(); \
} \
 \
/* some utility routines */  \
inline long deg(const TX& a) { \
  return a.rep.length()-1; \
} \
 \
/* Returns a read-only reference to coefficient i (zero if i not in range). \
 */ \
inline const T& coeff(const TX& a, long i) { \
  return a[i]; \
} \
 \
/* Read-only reference to leading term of a (zero if a==0). */ \
inline const T& LeadCoeff(const TX& a) { \
  return a[deg(a)]; \
} \
/* Read-only reference to constant term of a (zero if a==0). */ \
inline const T& ConstTerm(const TX& a) { \
  return a[0]; \
} \
 \
/* Makes coefficient of X^i equal to a. */ \
void SetCoeff(TX& x, long i, const T& a); \
 \
inline void GetCoeff(T& x, TX& a, long i) { \
  x=a[i]; \
} \
 \
void VectorCopy(vec_T& x, const TX& a, long n); \
 \
inline vec_T VectorCopy(const TX& a, long n) { \
  vec_T x; VectorCopy(x,a,n); NTL_OPT_RETURN(vec_T, x); \
} \
 \
/* x = reverse of a[0]..a[hi] (hi >= -1) */ \
void reverse(TX& x, const TX& a, long hi); \
inline TX reverse(const TX& a, long hi) { \
  TX x; reverse(x,a,hi); NTL_OPT_RETURN(TX, x); \
} \
 \
/* hi defaults to deg(a) */ \
inline void reverse(TX& x, const TX& a) { \
  reverse(x,a,deg(a)); \
} \
inline TX reverse(const TX& a) { \
  TX x; reverse(x,a,deg(a)); NTL_OPT_RETURN(TX, x); \
} \
 \
/* shift operators */   \
void LeftShift(TX& x, const TX& a, long n); \
 \
inline TX LeftShift(const TX& a, long n) { \
  TX x; LeftShift(x,a,n); NTL_OPT_RETURN(TX,x); \
} \
 \
void RightShift(TX& x, const TX& a, long n); \
 \
inline TX RightShift(const TX& a, long n) { \
  TX x; RightShift(x,a,n); NTL_OPT_RETURN(TX,x); \
} \
 \
inline TX& operator<<=(TX& x, long n) { \
  LeftShift(x,x,n); \
  return x; \
} \
 \
inline TX& operator>>=(TX& x, long n) { \
  RightShift(x,x,n); \
  return x; \
} \
 \
inline TX operator<<(const TX& a, long n) { \
  TX x; LeftShift(x,a,n); NTL_OPT_RETURN(TX, x); \
} \
 \
inline TX operator>>(const TX& a, long n) { \
  TX x; RightShift(x,a,n); NTL_OPT_RETURN(TX, x); \
} \
 \
/* arithmetic mod X^n */ \
/* x = a % X^m */ \
void trunc(TX& x, const TX& a, long m);  \
inline TX trunc(const TX& a, long m) { \
  TX x; trunc(x,a,m); NTL_OPT_RETURN(TX, x); \
} \



#define NTL_eq_polynomial_decl(T,vec_T,TX)  \
inline bool operator==(const TX& a, const TX& b) { \
  return a.rep==b.rep; \
} \
inline bool operator!=(const TX& a, const TX& b) { \
  return a.rep!=b.rep; \
} \
 \
/* promotion from T */ \
inline bool operator==(const T& a, const TX& b) { \
  return (deg(b)==0&&a==b[0]) || (deg(b)==-1&&IsZero(a)); \
} \
inline bool operator==(const TX& a, const T& b) { \
  return (deg(a)==0&&a[0]==b) || (deg(a)==-1&&IsZero(b)); \
} \
inline bool operator!=(const T& a, const TX& b) { \
  if (deg(b)==-1) return !IsZero(a); \
  return deg(b)>0 || a!=b[0]; \
} \
inline bool operator!=(const TX& a, const T& b) { \
  if (deg(a)==-1) return !IsZero(b); \
  return deg(a)>0 || a[0]!=b; \
} \



#define NTL_io_polynomial_decl(T,vec_T,TX)  \
NTL_SNS istream& operator>>(NTL_SNS istream&, TX&); \
NTL_SNS ostream& operator<<(NTL_SNS ostream&, const TX&); \




// addition, subtraction and negation methods
#define NTL_add_polynomial_decl(T,vec_T,TX)  \
void add(TX& x, const TX& a, const TX& b); \
void sub(TX& x, const TX& a, const TX& b); \
void negate(TX& x, const TX& a); \
 \
inline TX operator+(const TX& a, const TX& b) { \
  TX x; add(x,a,b); NTL_OPT_RETURN(TX, x); \
} \
inline TX operator-(const TX& a, const TX& b) { \
  TX x; sub(x,a,b); NTL_OPT_RETURN(TX, x); \
} \
inline TX operator-(const TX& a) { \
  TX x; negate(x,a); NTL_OPT_RETURN(TX, x); \
} \
 \
inline TX& operator+=(TX& x, const TX& a) { \
  add(x, x, a);  return x; \
} \
 \
inline TX& operator-=(TX& x, const TX& a) { \
  sub(x, x, a);  return x; \
} \
 \
/* promotions from T */ \
void add(TX& x, const T& a, const TX& b); \
void add(TX& x, const TX& a, const T& b); \
void sub(TX& x, const T& a, const TX& b); \
void sub(TX& x, const TX& a, const T& b); \
 \
inline TX operator+(const T& a, const TX& b) { \
  TX x; add(x,a,b); NTL_OPT_RETURN(TX,x); \
} \
inline TX operator+(const TX& a, const T& b) { \
  TX x; add(x,a,b); NTL_OPT_RETURN(TX,x); \
} \
inline TX operator-(const T& a, const TX& b) { \
  TX x; sub(x,a,b); NTL_OPT_RETURN(TX,x); \
} \
inline TX operator-(const TX& a, const T& b) { \
  TX x; sub(x,a,b); NTL_OPT_RETURN(TX,x); \
} \
 \
inline TX& operator+=(TX& x, const T& a) { \
  add(x,x,a); return x; \
} \
 \
inline TX& operator-=(TX& x, const T& a) { \
  sub(x,x,a); return x; \
} \


// scalar and polynomial multiplication
#define NTL_mul_polynomial_decl(T,vec_T,TX)  \
inline void set(TX& x) { \
  x.rep.SetLength(1); \
  set(x.rep.RawGet(0)); \
} \
inline bool IsOne(const TX& x) { \
  return x.rep.length()==1 && IsOne(x.rep.RawGet(0)); \
} \
inline void SetX(TX& x) { \
  x.rep.SetLength(2); \
  clear(x.rep.RawGet(0)); \
  set(x.rep.RawGet(1)); \
} \
inline bool IsX(const TX& x) { \
  return x.rep.length()==2 \
    && IsZero(x.rep.RawGet(0)) && IsOne(x.rep.RawGet(1)); \
} \
inline bool IsMonic(const TX& x) { \
  return x.rep.length()>0 && IsOne(x.rep.RawGet(x.rep.length()-1)); \
} \
 \
void SetCoeff(TX& x, long i); \
 \
void mul(TX& x, const TX& a, const T& b); \
void mul(TX& x, const T& a, const TX& b); \
 \
inline TX operator*(const T& a, const TX& b) { \
  TX x; mul(x,a,b); NTL_OPT_RETURN(TX, x); \
} \
inline TX operator*(const TX& a, const T& b) { \
  TX x; mul(x,a,b); NTL_OPT_RETURN(TX, x); \
} \
inline TX operator*=(TX& x, const T& a) { \
  mul(x,x,a); \
  return x; \
} \
 \
void mul(TX& x, const TX& a, const TX& b); \
void mul_classic(TX& x, const TX& a, const TX& b); \
inline TX operator*(const TX& a, const TX& b) { \
  TX x;  mul(x,a,b);  NTL_OPT_RETURN(TX,x); \
} \
inline TX operator*=(TX& x, const TX& a) { \
  mul(x,x,a); \
  return x; \
} \
void sqr(TX& x, const TX& a); \
void sqr_classic(TX& x, const TX& a); \
inline TX sqr(const TX& a) { \
  TX x;  sqr(x,a);  NTL_OPT_RETURN(TX,x); \
} \
 \
void MulTrunc(TX& x, const TX& a, const TX& b, long n); \
inline TX MulTrunc(const TX& a, const TX& b, long n) { \
  TX x;  MulTrunc(x,a,b,n);  NTL_OPT_RETURN(TX,x); \
} \
 \
void SqrTrunc(TX& x, const TX& a, long n); \
inline TX SqrTrunc(const TX& a, long n) { \
  TX x;  SqrTrunc(x,a,n);  NTL_OPT_RETURN(TX,x); \
} \
 \
void power(TX& x, const TX& a, long e); \
inline TX power(const TX& a, long e) { \
  TX x;  power(x,a,e);  NTL_OPT_RETURN(TX,x); \
} \
 \
void eval(T& x, const TX& f, const T& a); \
inline T eval(const TX& f, const T& a) { \
  T x; eval(x,f,a); return x; \
} \
 \
void BuildFromRoots(TX& x, const vec_T& a); \
inline TX BuildFromRoots(const vec_T& a) { \
  TX x;  BuildFromRoots(x,a);  NTL_OPT_RETURN(TX,x); \
} \




#define NTL_diff_polynomial_decl(T,vec_T,TX) \
void diff(TX& x, const TX& a); \
inline TX diff(const TX& a) { \
  TX x;  diff(x,a);  NTL_OPT_RETURN(TX,x); \
} \



#define NTL_char_polynomial_decl(T,vec_T,mat_T,TX) \
/* f = characteristic polynomial of M */ \
void CharPoly(TX& f, const mat_T& M); \



/**************** implementation macros ****************/


// zero_T is a variable containing the value 0 that can be 
// (read-only) referenced
#define NTL_polynomial_impl(T,vec_T,TX,zero_T) \
TX::TX(const TX& other, bool _normalize) {  \
  rep = other.rep; \
  if (_normalize) normalize(); \
} \
 \
void SetCoeff(TX& x, long i, const T& a) { \
  long len=x.rep.length(); \
  if (i>=len) { \
    if (!IsZero(a)) { \
      /* a might alias an element of x.rep so we need a copy */ \
      T aa(a); \
      x.ExtendTo(i+1); \
      /* swap is faster than assignment */ \
      swap(x.rep.RawGet(i),aa); \
    } \
  } \
  else if (i==len-1 && IsZero(a)) { \
    x.rep.SetLength(i); \
    while (i>0 && IsZero(x.rep.RawGet(--i))) \
      x.rep.SetLength(i); \
  } \
  else \
    x.rep.RawGet(i)=a; \
} \
 \
void VectorCopy(vec_T& x, const TX& a, long n) { \
  x.SetLength(n); \
  while (--n>=0) \
    x.RawGet(n) = a[n]; \
} \
 \
void reverse(TX& x, const TX& a, long hi) { \
  if (&x==&a) { \
    if (hi<x.rep.length()) \
      x.rep.SetLength(hi+1); \
    else \
      x.ExtendTo(hi+1); \
    for (long i=(hi-1)/2; i>=0; --i) \
      swap(x.rep.RawGet(i),x.rep.RawGet(hi-i)); \
  } \
  else { \
    x.rep.SetLength(hi+1); \
    for (long i=0; i<=hi; ++i) \
      x.rep.RawGet(i) = a[hi-i]; \
  } \
  x.normalize(); \
} \
 \
void LeftShift(TX& x, const TX& a, long n) { \
  if (n<=0) { \
    if (n<0) RightShift(x,a,-n); \
    else if (&x!=&a) x=a; \
    return; \
  } \
  long an = a.rep.length(); \
  if (&x==&a) { \
    x.ExtendTo(an+n); \
    for (long i=an-1; i>=0; --i) /* swap is faster than assignment */ \
      swap(x.rep.RawGet(i+n),x.rep.RawGet(i)); \
  } \
  else { \
    x.rep.SetLength(an+n); \
    for (long i=an-1; i>=0; --i) \
      x.rep.RawGet(i+n)=a.rep.RawGet(i); \
    do clear(x.rep.RawGet(--n)); while (n>=0); \
  } \
} \
 \
void RightShift(TX& x, const TX& a, long n) { \
  if (n<=0) { \
    if (n<0) LeftShift(x,a,-n); \
    else if (&x!=&a) x=a; \
    return; \
  } \
  if (n>=a.rep.length()) { \
    x.rep.SetLength(0); \
    return; \
  } \
  long xn = a.rep.length()-n; \
  if (&x==&a) { \
    for (long i=0; i<xn; ++i) /* swap is faster than assignment */ \
      swap(x.rep.RawGet(i),x.rep.RawGet(i+n)); \
    x.rep.SetLength(xn); \
  } \
  else { \
    x.rep.SetLength(xn); \
    for (long i=0; i<xn; ++i) \
      x.rep.RawGet(i)=a.rep.RawGet(i+n); \
  } \
} \
 \
void trunc(TX& x, const TX& a, long m) { \
  if (m<a.rep.length()) { \
    if (&x==&a) \
      x.rep.SetLength(m); \
    else \
      VectorCopy(x.rep,a,m); \
    x.normalize(); \
  } \
  else if (&x!=&a) \
    x.rep=a.rep; \
} \



#define NTL_eq_polynomial_impl(T,vec_T,TX) \
/* all methods inline, reserved for future use */ \



#define NTL_io_polynomial_impl(T,vec_T,TX) \
NTL_SNS istream & operator>>(NTL_SNS istream& in, TX& a) { \
  in>>a.rep; \
  a.normalize(); \
  return in; \
} \
  \
NTL_SNS ostream& operator<<(NTL_SNS ostream& out, const TX& a) { \
  out<<a.rep; \
  return out; \
} \



#define NTL_add_polynomial_impl(T,vec_T,TX)  \
void add(TX& x, const TX& a, const TX& b) { \
  if (&x==&b) { \
    /* x+=a */ \
    if (x.rep.length()<a.rep.length()) \
      x.ExtendTo(a.rep.length()); \
    for (long i=x.rep.length()-1; i>=0; --i) \
      add(x.rep.RawGet(i),a.rep.RawGet(i),b.rep.RawGet(i)); \
  } \
  else if (&x==&a) { \
    /* x+=b */ \
    if (x.rep.length()<b.rep.length()) \
      x.ExtendTo(b.rep.length()); \
    for (long i=x.rep.length()-1; i>=0; --i) \
      add(x.rep.RawGet(i),a.rep.RawGet(i),b.rep.RawGet(i)); \
  } \
  else { \
    long i = max(a.rep.length(),b.rep.length()); \
    x.rep.SetLength(i); \
    while (i>a.rep.length()) { \
      --i; \
      x.rep.RawGet(i) = b.rep.RawGet(i); \
    } \
    while (i>b.rep.length()) { \
      --i; \
      x.rep.RawGet(i) = a.rep.RawGet(i); \
    } \
    while (i>0) { \
      --i; \
      add(x.rep.RawGet(i),a.rep.RawGet(i),b.rep.RawGet(i)); \
    } \
  } \
  x.normalize(); \
} \
void sub(TX& x, const TX& a, const TX& b) { \
  if (&x==&b) { \
    /* x=a-x */ \
    if (x.rep.length()<a.rep.length()) \
      x.ExtendTo(a.rep.length()); \
    for (long i=x.rep.length()-1; i>=0; --i) \
      sub(x.rep.RawGet(i),a.rep.RawGet(i),b.rep.RawGet(i)); \
  } \
  else if (&x==&a) { \
    /* x-=b */ \
    if (x.rep.length()<b.rep.length()) \
      x.ExtendTo(b.rep.length()); \
    for (long i=x.rep.length()-1; i>=0; --i) \
      sub(x.rep.RawGet(i),a.rep.RawGet(i),b.rep.RawGet(i)); \
  } \
  else { \
    long i = max(a.rep.length(),b.rep.length()); \
    x.rep.SetLength(i); \
    while (i>a.rep.length()) { \
      --i; \
      negate(x.rep.RawGet(i),b.rep.RawGet(i)); \
    } \
    while (i>b.rep.length()) { \
      --i; \
      x.rep.RawGet(i) = a.rep.RawGet(i); \
    } \
    while (i>0) { \
      --i; \
      sub(x.rep.RawGet(i),a.rep.RawGet(i),b.rep.RawGet(i)); \
    } \
  } \
  x.normalize(); \
} \
void negate(TX& x, const TX& a) { \
  long an = a.rep.length(); \
  if (&x!=&a) x.rep.SetLength(an); \
  for (long i=0; i<an; ++i) \
    negate(x.rep.RawGet(i),a.rep.RawGet(i)); \
} \
/* promotions from T */ \
void add(TX& x, const T& a, const TX& b) { \
  long bn = b.rep.length(); \
  if (bn>1) { \
    if (&x!=&b) { \
      /* a could be an element of x */ \
      T sum; \
      add(sum,a,b.rep.RawGet(0)); \
      x=b; \
      swap(x.rep.RawGet(0),sum); \
    } \
    else \
      add(x.rep.RawGet(0),a,b.rep.RawGet(0)); \
  } \
  else if (bn==1) { \
    x.rep.SetLength(1); \
    add(x.rep.RawGet(0),a,b.rep.RawGet(0)); \
    if (IsZero(x.rep.RawGet(0))) clear(x); \
  } \
  else if (!IsZero(a)) \
    x=a; \
  else \
    clear(x); \
} \
void add(TX& x, const TX& a, const T& b) { \
  long an = a.rep.length(); \
  if (an>1) { \
    if (&x!=&a) { \
      /* b could be an element of x */ \
      T sum; \
      add(sum,a.rep.RawGet(0),b); \
      x=a; \
      swap(x.rep.RawGet(0),sum); \
    } \
    else \
      add(x.rep.RawGet(0),a.rep.RawGet(0),b); \
  } \
  else if (an==1) { \
    x.rep.SetLength(1); \
    add(x.rep.RawGet(0),a.rep.RawGet(0),b); \
    if (IsZero(x.rep.RawGet(0))) clear(x); \
  } \
  else if (!IsZero(b)) \
    x=b; \
  else \
    clear(x); \
} \
void sub(TX& x, const T& a, const TX& b) { \
  long bn = b.rep.length(); \
  if (bn>1) { \
    /* a could be an element of x */ \
    T diff; \
    sub(diff,a,b.rep.RawGet(0)); \
    negate(x,b); \
    swap(x.rep.RawGet(0),diff); \
  } \
  else if (bn==1) { \
    x.rep.SetLength(1); \
    sub(x.rep.RawGet(0),a,b.rep.RawGet(0)); \
    if (IsZero(x.rep.RawGet(0))) clear(x); \
  } \
  else if (!IsZero(a)) \
    x=a; \
  else \
    clear(x); \
} \
void sub(TX& x, const TX& a, const T& b) { \
  long an = a.rep.length(); \
  if (an>1) { \
    if (&x!=&a) { \
      /* b could be an element of x */ \
      T diff; \
      sub(diff,a.rep.RawGet(0),b); \
      x=a; \
      swap(x.rep.RawGet(0),diff); \
    } \
    else \
      sub(x.rep.RawGet(0),a.rep.RawGet(0),b); \
  } \
  else if (an==1) { \
    x.rep.SetLength(1); \
    sub(x.rep.RawGet(0),a.rep.RawGet(0),b); \
    if (IsZero(x.rep.RawGet(0))) clear(x); \
  } \
  else if (!IsZero(b)) { \
    x.rep.SetLength(1); \
    negate(x.rep.RawGet(0),b); \
  } \
  else \
    clear(x); \
} \




// scalar multiplication
#define NTL_mul_polynomial_impl_plain(T,vec_T,TX)  \
void SetCoeff(TX& x, long i) { \
  x.ExtendTo(i+1); \
  set(x.rep.RawGet(i)); \
} \
 \
void mul(TX& x, const T& a, const TX& b) { \
  if (IsZero(a)) \
    clear(x); \
  else { \
    T aa(a); /* a could be a coefficient in x */ \
    if (&x!=&b) x.rep.SetLength(b.rep.length()); \
    for (long i=0; i<b.rep.length(); ++i) \
      mul(x.rep.RawGet(i),aa,b.rep.RawGet(i)); \
    x.normalize(); \
  } \
} \
void mul(TX& x, const TX& a, const T& b) { \
  if (IsZero(b)) \
    clear(x); \
  else { \
    T bb(b); /* b could be a coefficient in x */ \
    if (&x!=&a) x.rep.SetLength(a.rep.length()); \
    for (long i=0; i<a.rep.length(); ++i) \
      mul(x.rep.RawGet(i),a.rep.RawGet(i),bb); \
    x.normalize(); \
  } \
} \
 \
void sqr_classic(TX& x, const TX& a) { \
  long da = deg(a); \
 \
  if (da<=0) { \
    if (da<0) \
      clear(x); \
    else { \
      x.rep.SetLength(1); \
      sqr(x.rep.RawGet(0),a.rep.RawGet(0)); \
    } \
    return; \
  } \
 \
  TX la; \
  const T *ap; \
  if (&x==&a) { \
    la = a; \
    ap = la.rep.elts(); \
  } \
  else \
    ap = a.rep.elts(); \
 \
  long d = 2*da; \
  x.rep.SetLength(d+1); \
  T *xp = x.rep.elts(); \
 \
  T t; \
  for (long i=0; i<=d; ++i) { \
    long jmin = max(0,i-da); \
    long jmax = min(da,i); \
    long m = jmax-jmin+1; \
    long m2 = m>>1; \
    jmax = jmin+m2-1; \
    clear(xp[i]); \
    for (long j=jmin; j<=jmax; ++j) { \
      mul(t,ap[j],ap[i-j]); \
      add(xp[i],xp[i],t); \
    } \
    add(xp[i],xp[i],xp[i]); \
    if (m&1) { \
      sqr(t,ap[jmax+1]); \
      add(xp[i],xp[i],t); \
    } \
  } \
 \
  x.normalize(); \
} \
 \
void mul_classic(TX& x, const TX& a, const TX& b) { \
  if (&a==&b) { \
    sqr_classic(x,a); \
    return; \
  } \
 \
  long da = deg(a); \
  long db = deg(b); \
 \
  if (da<0 || db<0) { \
    clear(x); \
    return; \
  } \
 \
  TX la, lb; \
  const T *ap, *bp; \
  if (&x==&a) { \
    la = a; \
    ap = la.rep.elts(); \
  } \
  else \
    ap = a.rep.elts(); \
 \
  if (&x==&b) { \
    lb = b; \
    bp = lb.rep.elts(); \
  } \
  else \
    bp = b.rep.elts(); \
 \
  long d = da+db; \
  x.rep.SetLength(d+1); \
  T *xp = x.rep.elts(); \
 \
  T t; \
  for (long i=0; i<=d; ++i) { \
    clear(xp[i]); \
    long jmax = min(da,i); \
    for (long j=max(0,i-db); j<=jmax; ++j) { \
      mul(t,ap[j],bp[i-j]); \
      add(xp[i],xp[i],t); \
    } \
  } \
  x.normalize(); \
} \
 \
void MulTrunc(TX& x, const TX& a, const TX& b, long n) { \
  TX t; \
  mul(t,a,b); \
  trunc(x,t,n); \
} \
 \
void SqrTrunc(TX& x, const TX& a, long n) { \
  TX t; \
  sqr(t,a); \
  trunc(x,t,n); \
} \
 \
void power(TX& x, const TX& a, long e) { \
  if (e<0) \
    Error("power: negative exponent"); \
 \
  if (e==0) { \
    set(x); \
    return; \
  } \
 \
  if (e==1 || IsZero(a) || IsOne(a)) { \
    if (&x!=&a) x=a; \
    return; \
  } \
 \
  long da = deg(a); \
 \
  /* Is power(T,long) available? \
  if (da==0) { \
    x = power(ConstTerm(a),e); \
    return; \
  } \
  */ \
 \
  if (da>(NTL_MAX_LONG-1)/e) \
    Error("overflow in power"); \
 \
  x.SetMaxLength(da*e+1); \
  set(x); \
 \
  TX ca(a); \
 \
  long k = NumBits(e); \
  for (long i=k-1; i>=0; --i) { \
    sqr(x,x); \
    if (bit(e,i)) \
      mul(x,x,ca); \
  } \
} \
 \
void eval(T& x, const TX& f, const T& a) { \
  /* Horner evaluation */ \
  T aa(a); \
  x = LeadCoeff(f); \
  for (long i=deg(f)-1; i>=0; --i) { \
    mul(x,x,aa); \
    add(x,x,f.rep.RawGet(i)); \
  } \
} \
 \
void IterBuild(T* a, long n) { \
  if (n<=0) return; \
  T b,t; \
  negate(a[0],a[0]); \
  for (long k=1; k<=n-1; ++k) { \
    negate(b,a[k]); \
    add(a[k],b,a[k-1]); \
    for (long i=k-1; i>=1; --i) { \
      mul(t,a[i],b); \
      add(a[i],t,a[i-1]); \
    } \
    mul(a[0],a[0],b); \
  } \
} \
void BuildFromRoots_classic(TX& x, const vec_T& a) { \
  long n = a.length(); \
  if (n==0) { \
    set(x); \
    return; \
  } \
 \
  x.rep.SetMaxLength(n+1); \
  x.rep = a; \
  IterBuild(x.rep.elts(),n); \
  x.rep.SetLength(n+1); \
  SetCoeff(x,n); \
} \



#define NTL_mul_polynomial_impl(T,vec_T,TX) \
NTL_mul_polynomial_impl_plain(T,vec_T,TX); \
void mul(TX& x, const TX& a, const TX& b) { \
  mul_classic(x,a,b); \
} \
void sqr(TX& x, const TX& a) { \
  sqr_classic(x,a); \
} \
void BuildFromRoots(TX& x, const vec_T& a) { \
  BuildFromRoots_classic(x,a); \
} \



#define NTL_diff_polynomial_impl(T,vec_T,TX) \
void diff(TX& x, const TX& a) { \
  long an = a.rep.length(); \
  if (an<=1) \
    clear(x); \
  else { \
    /* NOTE: if &x==&a, we are indexing beyond vector a, but since vector \
     * doesn't reallocate or destroy objects when being shrunk, it's ok. */ \
    x.rep.SetLength(an-1); \
    x.rep.RawGet(0) = a.rep.RawGet(1); \
    for (long i=2; i<an; ++i) \
      mul(x.rep.RawGet(i-1),i,a.rep.RawGet(i)); \
  } \
} \



#define NTL_char_polynomial_impl(T,vec_T,mat_T,TX) \
void CharPoly(TX& f, const mat_T& M) { \
  long n = M.NumRows(); \
  if (M.NumCols()!=n) \
    Error("CharPoly: nonsquare matrix"); \
 \
  if (n == 0) { \
    set(f); \
    return; \
  } \
 \
  T t; \
 \
  if (n == 1) { \
    SetX(f); \
    negate(t,M(1,1)); \
    SetCoeff(f,0,t); \
    return; \
  } \
 \
  mat_T H; \
  H = M; \
 \
  T u,t1; \
 \
  for (long m=2; m<=n-1; ++m) { \
    long i=m; \
    while (i<=n && IsZero(H(i,m-1))) \
      ++i; \
    if (i<=n) { \
      t = H(i,m-1); \
      if (i>m) { \
        swap(H(i),H(m)); \
        /* swap columns i and m */ \
        for (long j=1; j<=n; ++j) \
          swap(H(j,i),H(j,m)); \
      } \
      for (i=m+1; i<=n; ++i) { \
        div(u,H(i,m-1),t); \
        for (long j=m; j<=n; ++j) { \
          mul(t1,u,H(m,j)); \
          sub(H(i,j),H(i,j),t1); \
        } \
        for (long j=1; j<=n; ++j) { \
          mul(t1,u,H(j,i)); \
          add(H(j,m),H(j,m),t1); \
        } \
      } \
    } \
  } \
 \
  TX F[n+1]; \
  TX T1; \
  T1.SetMaxLength(n); \
 \
  set(F[0]); \
  for (long m=1; m<=n; ++m) { \
    LeftShift(F[m],F[m-1],1); \
    mul(T1,F[m-1],H(m,m)); \
    sub(F[m],F[m],T1); \
 \
    set(t); \
    for (long i=1; i<=m-1; ++i) { \
      mul(t,t,H(m-i+1,m-i)); \
      mul(t1,t,H(m-i,m)); \
      mul(T1,F[m-i-1],t1); \
      sub(F[m],F[m],T1); \
    } \
  } \
  f = F[n]; \
} \


#endif
