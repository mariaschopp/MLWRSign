#include <cstring>
#include <cstdlib>

#include "verify.h"
#include "sign.h"
#include "keygen.h"
#include <cstdio>

bool verify (vk_t vk, unsigned char *m, unsigned int mlen, signat_t sig) {
    unsigned char mu[CRHBYTES];
    polyveck w1_prime;
    poly c_prime;

    // Calculate mu
    calculate_mu (mu, vk, m, mlen);

    // Calculate w' from hint vector and signature
    calculate_w1_prime (w1_prime, vk, sig);

    // Calculate c' = H(mu, w1')
    memset(&c_prime, 0, sizeof(poly));
    challenge (c_prime, mu, w1_prime);

    // Check acceptance criteria.
    return verify_signature(sig, c_prime);
}

int verify_from_bytes (unsigned char *m, unsigned long long *mlen,
                       unsigned char *sm, unsigned long long smlen, unsigned char * pk)
{
    vk_t vk;
    signat_t sig;

    #ifdef MATRIX_A
    vk.A = matrix_A;
    #else
    fprintf (stderr, "MATRIX_A undefined. Public A must be defined for verifying from key bytes.\n");
    return -1;
    #endif

    deserialize_verifykey (vk, pk);
    *mlen = unpack_signed_message (sm, sig, m, smlen);

    if (verify (vk, m, (unsigned int) *mlen, sig))
        return 0;

    return -3;
}

// SUBFUNCTIONS

void calculate_mu (unsigned char *mu, vk_t vk, unsigned char *m, unsigned int mlen) {
    unsigned char tr[CRHBYTES];
    int buflen = mlen + CRHBYTES;
    unsigned char buf[buflen];

    calculate_tr (tr, vk.t1);

    // mu = hash (tr || m)
    memcpy (buf, tr, CRHBYTES);
    memcpy (buf + CRHBYTES, m, mlen);
    shake256 (mu, CRHBYTES, buf, mlen + CRHBYTES);
}

void calculate_w1_prime (polyveck &w1_prime, vk_t vk, signat_t sig) {
    polyveck az, ct1, r;
    uint32_t i, j;

    // Calculate az = round(A * z)
    mult_MatVecl(az, vk.A, sig.z, PQS_q);

    for (i=0; i<PQS_k; i++) {
        for (j=0; j<PQS_n; j++) {
            az.polynomial[i].coeffs[j] = round (az.polynomial[i].coeffs[j]);
        }
    }

    // Calculate ct1 = c * t1 * 2^d
    for (i = 0; i < PQS_k; i++) {
        MULT_POLYPOLY (sig.c, vk.t1.polynomial[i], ct1.polynomial[i], PQS_n, PQS_q);

        for (j = 0; j < PQS_n; j++) {
            ct1.polynomial[i].coeffs[j] <<= PQS_d;
        }
    }

    // Calculate r = az - ct1
    for (i = 0; i < PQS_k; i++) {
        subt_PolyPoly (r.polynomial[i], az.polynomial[i], ct1.polynomial[i], PQS_q);
    }

    use_hints (w1_prime, sig.h, r, PQS_gamma_bar2_bits + 1);
}

void use_hints (polyveck &w1_prime, polyveck h, polyveck r, uint32_t alphabits) {
    uint32_t i, j, m, r1ij;
    polyveck r1, r0;

    m = 1 << (PQS_p_bits - alphabits);
    decompose (r1, r0, r, alphabits);

    for (i = 0; i < PQS_k; i++) {
        for (j = 0; j < PQS_n; j++) {
            r1ij = r1.polynomial[i].coeffs[j];
            w1_prime.polynomial[i].coeffs[j] = r1ij; 

            if (h.polynomial[i].coeffs[j] == 1) {
                if (r0.polynomial[i].coeffs[j] > 0) {
                    w1_prime.polynomial[i].coeffs[j] = mod_neg (r1ij + 1, m);
                }
                else {
                    w1_prime.polynomial[i].coeffs[j] = mod_neg (r1ij - 1, m);
                }
            }
        }
    }
}

bool verify_signature(signat_t sig, poly c_prime) {
    uint32_t i, j, bound, hw;

    // Verify criterion 1: z inf-norm bound
    bound = PQS_gamma1 - PQS_beta1;

    for (i=0; i<PQS_l; i++) {
        for (j=0; j<PQS_n; j++){
            if (abs(sig.z.polynomial[i].coeffs[j]) > bound){
                return false;
            }
        }
    }

    // Verify criterion 2: c = c_prime
    for (j=0; j<PQS_n ; j++) {
        if (sig.c.coeffs[j] != c_prime.coeffs[j]) {
            return false;
        }
    }

    // Verify criterion 3: hw(h) <= omega
    hw = 0;
    for (i = 0; i < PQS_k; i++) {
        for (j = 0; j < PQS_n; j++) {
            if (sig.h.polynomial[i].coeffs[j] == 1) {
                hw++;
            }
        }
    }

    if (hw > PQS_omega) {
        return false;
    }

    return true;
}
