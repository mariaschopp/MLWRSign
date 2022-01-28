#ifndef VERIFY
#define VERIFY

#include "utils/utils.h"

bool verify(vk_t vk, unsigned char *m, unsigned int mlen, signat_t sig);

int verify_from_bytes (unsigned char *m, unsigned long long *mlen, 
                       unsigned char *sm, unsigned long long smlen, unsigned char * pk);

void calculate_mu (unsigned char *mu, vk_t vk, unsigned char *m, unsigned int mlen);
void calculate_w1_prime (polyveck &w1_prime, vk_t vk, signat_t sig);
void use_hints (polyveck &w1_prime, polyveck h, polyveck r, uint32_t alphabits);
bool verify_signature(signat_t sig, poly c_prime);

#endif