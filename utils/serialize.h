#ifndef PQS_J_SERIALIZE_H
#define PQS_J_SERIALIZE_H

#include "params.h"
#include "serialize_common.h"

void serialize_secret (sk_t &sk, unsigned char *sk_bytes);
void deserialize_secret (sk_t &sk, unsigned char *sk_bytes);

void serialize_verifykey (vk_t &vk, unsigned char *vk_bytes);
void deserialize_verifykey (vk_t &vk, unsigned char *vk_bytes);

unsigned long long pack_signed_message (unsigned char *sm, signat_t &sig, unsigned char *m, unsigned long long mlen);
unsigned long long unpack_signed_message (unsigned char *sm, signat_t &sig, unsigned char *m, unsigned long long smlen);

int sig_compare (signat_t sig1, signat_t sig2);
int poly_compare (poly p1, poly p2, int length);

#endif
