#ifndef ARITH_H
#define ARITH_H

#include "params.h"
#include <cstdint>

void mult_MatVecl(polyveck &res, polymatkl A, polyvecl v, int32_t modulus);
void add_PolyPoly(poly &res, poly p1, poly p2, int32_t modulus);
void subt_PolyPoly(poly &res, poly p1, poly p2, int32_t modulus);

int32_t round(int32_t x);
int32_t mod(int32_t x, int32_t n);
int32_t mod_neg(int32_t x, int32_t n);

#endif
