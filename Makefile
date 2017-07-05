CC		= gcc
CPP		= g++
LINK		= g++

NTLPREFIX	= /usr/local

PREFFLAGS	=
CPPFLAGS	= -I$(NTLPREFIX)/include -I. -O2 -Wall $(PREFFLAGS)
LINKFLAGS	= -L$(NTLPREFIX)/lib -lntl -lgmp $(PREFFLAGS)

ALL_PROGS	= ecpp-test
COMMON_OBJS	= EC_p.o ZZFactoring.o \
		  RRX.o CC.o HilbertClassPolynomials.o
COMMON_HEADERS	=
ALL_DIRS	=


all:	$(ALL_DIRS) $(ALL_PROGS)

complext.o:	/usr/include/g++-3/std/complext.cc
	$(CPP) $(CPPFLAGS) -c $^

ecpp-test:	ecpp-test.o $(COMMON_OBJS)
	$(LINK) -o $@ $^ $(LINKFLAGS)

clean:
	rm -f *% *~ *.o core a.out $(ALL_PROGS)
