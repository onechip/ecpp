
#include "RRX.h"
#include <NTL/new.h>


NTL_START_IMPL;

const RR RR_zero;

NTL_vector_impl(RR,vec_RR);
NTL_eq_vector_impl(RR,vec_RR);
NTL_io_vector_impl(RR,vec_RR);

NTL_polynomial_impl(RR,vec_RR,RRX,RR_zero);
NTL_eq_polynomial_impl(RR,vec_RR,RRX);
NTL_io_polynomial_impl(RR,vec_RR,RRX);
NTL_add_polynomial_impl(RR,vec_RR,RRX);
NTL_mul_polynomial_impl(RR,vec_RR,RRX);
NTL_diff_polynomial_impl(RR,vec_RR,RRX);



NTL_END_IMPL;
