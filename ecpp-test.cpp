#include <stdio.h>
#include <string.h>

#include <NTL/vec_vec_ZZ.h>
#include <NTL/ZZ_pXFactoring.h>
#include "ZZFactoring.h"
#include "HilbertClassPolynomials.h"
#include "EC_p.h"

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


int discriminant_array[] = { 
  // class number 1
  3,4,7,8,11,19,43,67,163,
  // class number 2
  15,20,24,35,40,51,52,88,91,115,123,148,187,232,235,267,403,427,
  // class number 3
  23,31,59,83,107,139,211,283,307,331,379,499,547,643,883,907,
  // class number 4
  //14,17,21,30,33,34,39,42,46,55,57,70,73,78,82,85,93,97,102,130,142,155,177,
  //190,193,195,203,219,253,259,291,323,355,435,483,555,595,627,667,715,723,
  //763,795,955,1003,1227,1243,1387,1411,1435,1507,1555, // and 133?
  // class number 5
  47,79,103,127,131,179,227,347,443,523,571,619,683,691,739,787,947,1051,1123,
  1723,1747,1867,2203,2347,2683,
  // class number 6
  //87,104,116,152,212,244,247,339,411,424,436,451,472,515,628,707,771,808,835,
  //843,856,1048,1059,1099,1108,1147,1192,1203,1219,1267,1315,1347,1363,1432,
  //1563,1588,1603,1843,1915,1963,2227,2283,2443,2515,2563,2787,2923,3235,3427,
  //3523,3763,
  // class number 7
  71,151,223,251,463,467,487,587,811,827,859,1163,1171,1483,1523,1627,1787,
  1987,2011,2083,2179,2251,2467,2707,3019,3067,3187,3907,4603,5107,5923,
  // class number 8
  // class number 9
  199,367,419,491,563,823,1087,1187,1291,1423,1579,2003,2803,3163,3259,3307,
  3547,3643,4027,4243,4363,4483,4723,4987,5443,6043,6427,6763,6883,7723,8563,
  8803,9067,10627,
  // class number 10
  // class number 11
  167,271,659,967,1283,1303,1307,1459,1531,1699,2027,2267,2539,2731,2851,2971,
  3203,3347,3499,3739,3931,4051,5179,5683,6163,6547,7027,7507,7603,7867,8443,
  9283,9403,9643,9787,10987,13003,13267,14107,14683,15667,
  // class number 12
  // class number 13
  191,263,607,631,727,1019,1451,1499,1667,1907,2131,2143,2371,2659,2963,3083,
  3691,4003,4507,4643,5347,5419,5779,6619,7243,7963,9547,9739,11467,11587,
  11827,11923,12043,14347,15787,16963,20563,
  // class number 14
  // class number 15
  239,439,751,971,1259,1327,1427,1567,1619,2243,2647,2699,2843,3331,3571,3803,
  4099,4219,5003,5227,5323,5563,5827,5987,6067,6091,6211,6571,7219,7459,7547,
  8467,8707,8779,9043,9907,10243,10267,10459,10651,10723,11083,11971,12163,
  12763,13147,13963,14323,14827,14851,15187,15643,15907,16603,16843,17467,
  17923,18043,18523,19387,19867,20707,22003,26203,27883,29947,32323,34483,
  // class number 16
  // class number 17
  383,991,1091,1571,1663,1783,2531,3323,3947,4339,4447,4547,4651,5483,6203,
  6379,6451,6827,6907,7883,8539,8731,9883,11251,11443,12907,13627,14083,14779,
  16699,17827,18307,19963,21067,23563,24907,25243,26083,26107,27763,31627,
  33427,36523,   //371213,
  // class number 18
  // class number 19
  311,359,919,1063,1543,1831,2099,2339,2459,3343,3463,3467,3607,4019,4139,4327,
  5059,5147,5527,5659,6803,8419,8923,8971,9619,10891,11299,15091,15331,16363,
  16747,17011,24763,26227,27043,29803,31123,37507,38707,
  // class number 20
  // class number 21
  431,503,743,863,1931,2503,2579,2767,2819,3011,3371,4283,4523,4691,5011,5647,
  5851,5867,6323,6691,7907,8059,8123,8171,8243,8387,8627,8747,9091,9187,9811,
  9859,10067,10771,11731,12107,12547,13171,13291,13339,13723,14419,14563,15427,
  16339,16987,17107,17707,17971,18427,18979,19483,19531,19819,20947,21379,
  22027,22483,22963,23227,23827,25603,26683,27427,28387,28723,28867,31963,
  32803,34147,34963,35323,36067,36187,39043,40483,44683,46027,49603,51283,
  52627,55603,58963,59467,61483,
  // class number 22
  // class number 23
  647,1039,1103,1279,1447,1471,1811,1979,2411,2671,3491,3539,3847,3923,4211,
  4783,5387,5507,5531,6563,6659,6703,7043,9587,9931,10867,10883,12203,12739,
  13099,13187,15307,15451,16267,17203,17851,18379,20323,20443,20899,21019,
  21163,22171,22531,24043,25147,25579,25939,26251,26947,27283,28843,30187,
  31147,31267,32467,34843,35107,37003,40627,40867,41203,42667,43003,45427,
  45523,47947,90787,
  // others ??
  12,16,27,28,32,36,39,44,48,55,56,60,63,64,68,72,75,76,80,84,87,92,95,96,99,
  // end of list
  0
};


/* Modified Cornacchia.  (Algorithm 1.5.3 in Cohen)
 * Requires: p prime, D negative and congruent to 0 or 1 (mod 4), and |D|<4p.
 * Returns true if x^2 + |D|y^2 = 4p, and false if no solution exists.
 */
bool Cornacchia(ZZ& x, ZZ& y, const ZZ& D, const ZZ& p) {
  if (p==2) {
    // step 1
    SqrRoot(x,D+8);
    set(y);
    return sqr(x)==D+8;
  }
  // step 2
  if (Jacobi(D,p)==-1) 
    return false;
  // step 3
  ZZ x0,Dm;
  rem(Dm,D,p);
  SqrRootMod(x0,Dm,p);
  if (IsOdd(D)!=IsOdd(x0))
    sub(x0,p,x0);
  ZZ a,b,l;
  a=2*p;
  b=x0;
  SqrRoot(l,4*p);
  // step 4
  while (b>l) {
    rem(a,a,b);
    swap(a,b);
  }
  // step 5
  ZZ c;
  if (!divide(c,4*p-sqr(b),-D))
    return false;
  ZZ cs;
  SqrRoot(cs,c);
  if (c!=sqr(cs))
    return false;
  x=b;
  y=cs;
  return true;
}



// find curve parameters
void find_curve(ZZ_p& a, ZZ_p& b, long D, const ZZ& N, 
		HilbertClassPolynomials& HCP) {
  if (D==-3) {
    clear(a);
    conv(b,-1);
    return;
  }
  else if (D==-4) {
    conv(a,-1);
    clear(b);
    return;
  }
  
  const ZZX& T = HCP.poly(D);
  cout<<"find_curve: T="<<T<<endl;
  ZZ_pX Tp;
  conv(Tp,T);
  cout<<"find_curve: Tp="<<Tp<<endl;
  ZZ_p j;
  FindRoot(j,Tp);
  cout<<"j="<<j<<endl;
  ZZ_p c;
  div(c,j,j-1728);
  mul(a,-3,c);
  mul(b,2,c);
}


// q is largest prime factor of m
// return true if q>=threshold and prime
bool check_for_factor(ZZ& q, const ZZ& m, const ZZ& threshold) {
  static const long bound=100000;

  // m must be composite
  if (ProbPrime(m))
    return false;

  // remove small factors
  q=m;
  ZZ t;
  PrimeSeq s;
  for (long p=s.next(); p<bound; p=s.next())
    while (DivRem(t,q,p)==0)
      q=t;

  if (q<threshold)
    return false;

  // lower bound on size of smallest factor
  // we are assuming that ECM finds factors smallest first (more or less)
  ZZ small_factor;
  conv(small_factor,bound);

  do {
    if (ProbPrime(q))
      return true;

    // q is composite and has a factor greater than small_factor,
    // thus it must be at least threshold*small_factor to have a 
    // large factor >threshold
    if (q<threshold*small_factor)
      return false;

    // attempt to factor with ECM
    // we are looking for factor such that
    //    small_factor < factor < q/threshold
    ECM(t,q,to_ZZ(NTL_MAX_LONG),0.5,true);

    if (IsOne(t))
      return false;

    if (t<threshold) {
      // new lower bound on remaining factors
      if (small_factor<t)
	small_factor=t;
      q/=t;
    }
    else {
      // this will be rare, but we should update small_factor anyway
      q=t;
    }

  } while (true);
}


/* 
 * Certificate format:
 *   [N q m/q a b x y]
 */
bool ProvePrime_Atkin(const ZZ& N, vec_ZZ& cert,
		      HilbertClassPolynomials& HCP) {
  ZZ threshold;
  // want: threshold = ceil((N^{1/2}+1)^2) = ceil(N^{1/2} + 2*N^{1/4} + 1)
  ZZ Ns,Nf;
  SqrRoot(Ns,N);    // Ns+1  = ceil(N^{1/2})
  SqrRoot(Nf,Ns+1); // Nf+1 >= ceil(N^{1/4})
  threshold = Ns+1 + 2*(Nf+1) + 1;  // threshold may be a little too high here

  cout<<endl;
  cout<<"N="<<N<<endl;
  cout<<"threshold = "<<threshold<<endl;

  ZZ_pBak bak1;
  bak1.save();
  ZZ_p::init(N);

  EC_pBak bak2;
  bak2.save();

  long n=0;
  ZZ D;
  while (!IsZero(D=-discriminant_array[n++])) {
    //    long D4 = rem(D,4);
    //if (D4!=0 && D4!=1)
    // continue;
    if (Jacobi(D+N,N)!=1)
      continue;
    // step 3
    ZZ u,v;
    if (!Cornacchia(u,v,D,N))
      continue;
    // step 4
    ZZ m,q;
    bool found=false;
    if (check_for_factor(q,m=N+1+u,threshold))
      found=true;
    else if (check_for_factor(q,m=N+1-u,threshold))
      found=true;
    else if (D==-4) {
      if (check_for_factor(q,m=N+1+2*v,threshold))
	found=true;
      else if (check_for_factor(q,m=N+1-2*v,threshold))
	found=true;
    }
    else if (D==-3) {
      if (check_for_factor(q,m=N+1+(u+3*v)/2,threshold))
	found=true;
      else if (check_for_factor(q,m=N+1-(u+3*v)/2,threshold))
	found=true;
      else if (check_for_factor(q,m=N+1+(u-3*v)/2,threshold))
	found=true;
      else if (check_for_factor(q,m=N+1-(u-3*v)/2,threshold))
	found=true;
    }
    if (!found)
      continue;
    cout<<"m="<<m<<endl;
    cout<<"q="<<q<<endl;
    // step 6
    EC_pCurve curve;
    ZZ_p A,B;
    find_curve(A,B,to_long(D),N,HCP);
    SetCoeff(curve,1,A);
    SetCoeff(curve,0,B);
    cout<<"curve="<<curve<<endl;
    // step 7
    ZZ_p g;
    do {
      random(g);
      if (Jacobi(rep(g),N)!=-1)
	continue;
      if ((N%3)==1) {
	// always happens in the case D==-3
	ZZ_p t;
	power(t,g,(N-1)/3);
	if (IsOne(t))
	  continue;
	//cout<<"t="<<t<<endl;
      }
      break;
    } while (true);
    //cout<<"g="<<g<<endl;
    long k=0;
    EC_p::init(curve);
    EC_p P,P1,P2;
    do {
      // step 8
      random(P);
      if (!IsValid(P)) 
	return false;
      // step 9
      //cout<<"P="<<P<<endl;
      mul(P2,P,m/q);
      //cout<<"P2="<<P2<<endl;
      // NOTE: if IsZero(P2) just choose a new point, not a new curve!!!
      if (!IsZero(P2)) {
	mul(P1,P2,q);
	//cout<<"P1="<<P1<<endl;
	if (IsZero(P1)) 
	  break;
      }
      // step 10
      ++k;
      if (D==-3) {
	if (k>=6) {
	  cout<<"P="<<P<<endl;
	  cout<<"P1="<<P1<<endl;
	  cout<<"P2="<<P2<<endl;
	  return false;
	}
	SetCoeff(curve,0,coeff(curve,0)*g);
      }
      else if (D==-4) {
	if (k>=4) {
	  cout<<"P="<<P<<endl;
	  cout<<"P1="<<P1<<endl;
	  cout<<"P2="<<P2<<endl;
	  return false;
	}
	SetCoeff(curve,1,coeff(curve,1)*g);
      }
      else {
	if (k>=2) {
	  cout<<"P="<<P<<endl;
	  cout<<"P1="<<P1<<endl;
	  cout<<"P2="<<P2<<endl;
	  return false;
	}
	SetCoeff(curve,1,coeff(curve,1)*sqr(g));
	SetCoeff(curve,0,coeff(curve,0)*sqr(g)*g);
      }
      EC_p::init(curve);
    } while (true);

    // check that P2 has an affine representation
    ZZ G;
    GCD(G,rep(P2.Z),N);
    if (!IsOne(G))
      return false;

    // step 13
    ZZ_p Px,Py;
    P.affine(Px,Py);
    cert.SetLength(7);
    cert[0] = N;
    cert[1] = q;
    cert[2] = m/q;
    cert[3] = coeff(curve,1)==-1 ? to_ZZ(-1) : rep(coeff(curve,1));
    cert[4] = coeff(curve,0)==-1 ? to_ZZ(-1) : rep(coeff(curve,0));
    cert[5] = rep(Px);
    cert[6] = rep(Py);
    return true;
  }
  Error("ProvePrime: ran out of discriminants");
  return false;
}


bool ProvePrime(const ZZ& N, vec_vec_ZZ& certs) {
  certs.SetLength(0);
  long i=0;
  ZZ num(N);
  HCP_generate HCP;
  while (NumBits(num)>30) {
    certs.SetLength(i+1);
    if (!ProvePrime_Atkin(num,certs[i],HCP))
      return false;
    num=certs[i][1];
    ++i;
  }
  // prove small prime by trial division
  long n,limit;
  conv(limit,sqrt(n));
  conv(n,num);
  PrimeSeq s;
  for (long p=s.next(); p<=limit; p=s.next())
    if (n%p==0)
      return false;
  return true;
}



void check_table() {
  int i=0;
  int max=0;
  char* used = new char[100000];
  memset(used,0,100000);
  while (discriminant_array[i]) {
    int d=discriminant_array[i];
    used[d]=1;
    if (d>max) max=d;
    if (d<3) 
      cerr<<"d="<<d<<"  i="<<i<<endl;
    long d4 = (-d)%4;
    if (d4<0) d4+=4;
    if (d4!=0 && d4!=1)
      cerr<<"d="<<d<<"  i="<<i<<"  d4="<<d4<<endl;
    int j=i+1;
    while (discriminant_array[j]) {
      if (d==discriminant_array[j])
	cerr<<"duplicate: d="<<d<<"  i="<<i<<"  j="<<j<<endl;
      ++j;
    }
    ++i;
  }
  cerr<<"total="<<i<<endl;
  cerr<<"max="<<max<<endl;
  for (int d=3; d<500; ++d)
    if (used[d]==0) {
      long d4 = (-d)%4;
      if (d4<0) d4+=4;
      if (d4==0 || d4==1)
	cout<<d<<",";
    }
  cout<<endl;
}

int main(int argc, char*argv[]) {
  //check_table();
  
  long bits=0;
  if (argc>1)
    bits=atoi(argv[1]);
  if (bits<32)
    bits=48;

  ZZ N;
  if (bits==89)
    N = (to_ZZ(1)<<89) - 1;
  else
    RandomPrime(N,bits);

  vec_vec_ZZ certs;
  if (!ProvePrime(N,certs)) {
    cout<<"Failed to prove prime.  Maybe N is composite?"<<endl;
    return 1;
  }

  cout<<endl;
  cout<<"N="<<N<<endl;
  for (long i=0; i<certs.length(); ++i) {
    cout<<endl;
    cout<<"[q="<<certs[i][1]<<endl;
    cout<<" r="<<certs[i][2]<<endl;
    cout<<" a="<<certs[i][3]<<endl;
    cout<<" b="<<certs[i][4]<<"]"<<endl;
  }
  
  //cout<<endl;
  //HCP_generate HCP;
  //cout<<"T(-499) = "<<HCP.poly(-499)<<endl;

  return 0;
}
