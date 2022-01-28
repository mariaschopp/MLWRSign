#include "arith.h"
#include "poly_mul.h"
#include <cstring>

// ARITHMETIC FUNCTIONS

// Multiply a k*l matrix by an l-vector, returning a k-vector.
void mult_MatVecl(polyveck &res, polymatkl A, polyvecl v, int32_t modulus) {
    int i, j;
    poly acc;
    memset(&res, 0, sizeof(polyveck));
    for (i=0; i<PQS_k; i++){
        for (j=0; j<PQS_l; j++){
            MULT_POLYPOLY(A.elem[i][j], v.polynomial[j], acc, PQS_n, modulus);
            add_PolyPoly (res.polynomial[i], res.polynomial[i], acc, modulus);
        }
    }
}

// Add polynomial p1 to polynomial p2.
void add_PolyPoly(poly &res, poly p1, poly p2, int32_t modulus) {
    for (int n=0; n<PQS_n; n++)
        res.coeffs[n] = mod_neg(p1.coeffs[n] + p2.coeffs[n], modulus);
}

// Subtract polynomial p2 from polynomial p1.
void subt_PolyPoly(poly &res, poly p1, poly p2, int32_t modulus) {
    for (uint32_t i=0; i<PQS_n; i++)
        res.coeffs[i] = mod_neg(p1.coeffs[i] - p2.coeffs[i], modulus);
}

// Round an integer PQS_q -> PQS_p
int32_t round(int32_t x) {
    int32_t h = 1 << (PQS_qdivp_bits - 1);
    int32_t sign_bit = 1 - 2*(x < 0);
    return sign_bit * ((x + h) >> (PQS_qdivp_bits));
}

// Reduce x to the range [0,n].
int32_t mod(int32_t x, int32_t n) {
    return x & (n - 1);
}

// Reduce x to the range (-2^n / 2, 2^n / 2].
int32_t mod_neg(int32_t x, int32_t n) {
    int32_t mn = 1 << n;
    int32_t y = mod(x, mn);
    return y - (y >> (n-1)) * mn;
}
