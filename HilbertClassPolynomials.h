#ifndef NTL_HilbertClassPolynomials__H
#define NTL_HilbertClassPolynomials__H

#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>

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


// abstract class for Hilbert class polynomials
class HilbertClassPolynomials {
private:
  vec_ZZX cache;

protected:
  // return reference to polynomial associated with D in cache
  // the polynomial will be zero if it hasn't been set yet
  ZZX& cached_poly(long D);

  /* Get or compute polynomial associated with negative discriminant D
   * and store in dest.  D%4 = 0 or 1.
   */
  virtual void get_poly(ZZX& dest, long D)=0;

public:
  virtual ~HilbertClassPolynomials() {}

  // get polynomial associated with negative discriminant D
  // if D>0, -D is used
  const ZZX& poly(long D);

  inline void poly_reduced(ZZ_pX& P, long D) {
    conv(P,poly(D));
  }
  
};


// class that computes polynomials
class HCP_generate : public HilbertClassPolynomials {
protected:
  void get_poly(ZZX& dest, long D);
};



#endif
