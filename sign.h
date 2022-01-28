#ifndef SIGN
#define SIGN

#include "utils/utils.h"

void sign (sk_t &sk, vk_t &vk, unsigned char *m, unsigned int mlen, signat_t &sig);
int sign_from_keybytes (unsigned char *sm, unsigned long long* smlen,
                        unsigned char *m, unsigned long long mlen,
                        unsigned char *sk, unsigned char *pk);

void gen_seeds (unsigned char *mu, unsigned char *seed, sk_t &sk, unsigned char *m, unsigned int mlen);
void calculate_s2 (polyveck &s2, polyveck t, polymatkl A, polyvecl s1);
void generate_y (polyvecl &y, unsigned char *seed);
void calculate_w_xi1 (polyveck &w, polyveck &xi1, polymatkl A, polyvecl y);
void calculate_w1 (polyveck &w1, polyveck w);
void calculate_z (polyvecl &z, polyvecl y, polyvecl s, poly c);
int validate_z_ct0 (polyvecl z, polyveck ct0, poly c, polyveck t0);
void calculate_xi2_nu (polyveck &nu, polyveck &xi2, polyveck &rcs2, poly c, polyveck s2, polyveck xi1);
void calculate_r (polyveck &r, polyveck w, polyveck rcs2, polyveck nu);
int validate_r (polyveck r0, polyveck r1, polyveck w1);
int make_hints (polyveck &h, polyveck ct0, polyveck r);

#endif