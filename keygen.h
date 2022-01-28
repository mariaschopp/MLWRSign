#ifndef KEYGEN_H
#define KEYGEN_H

#include "utils/utils.h"

void keygen(sk_t &sk, vk_t &vk);
int keygen_bytes (unsigned char *pk, unsigned char *sk);

void generate_A (polymatkl &A);
void generate_s1 (polyvecl &s);
void calculate_t (polyveck &t, polymatkl A, polyvecl s);
void calculate_tr (unsigned char *tr, polyveck t1);

#endif