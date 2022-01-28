#include <cstring>
#include <cstdlib>

#include "sign.h"
#include <cstdio>

void sign (sk_t &sk, vk_t &vk, unsigned char *m, unsigned int mlen, signat_t &sig) {
    // Note: s2, xi1, and xi2 would use floating point values which would undo the
    // performance benefit of choosing p and q to be powers of 2. Instead these are
    // calculated with an extra multiple of (q/p) to avoid using floats and the extra
    // multiples are accounted for later in the computation.

    polyveck t, s2, w, w1, xi1, xi2, rcs2, r, r1, r0, ct0, nu;
    polyvecl y;
    int hw;

    unsigned char seed[CRHBYTES];
    unsigned char mu[CRHBYTES];

    gen_seeds (mu, seed, sk, m, mlen);

    // Reconstruct t
    recompose (t, vk.t1, sk.t0, PQS_d);

    // Calculate 2^d*s2 = (2^d*t - A*s1)
    calculate_s2 (s2, t, vk.A, sk.s1);

    do {
        do {
            do {
                // Generate y from distribution
                generate_y (y, seed);

                // Calculate w = Round ((p/q) A*y), (q/p)*xi1 = (q/p)*w - A*y
                calculate_w_xi1 (w, xi1, vk.A, y);

                // Calculate w1 = HighBits(w)
                calculate_w1 (w1, w);

                // Calculate c = H(mu, w1)
                challenge (sig.c, mu, w1);

                // Calculate z = y + c*s1
                calculate_z (sig.z, y, sk.s1, sig.c);
            } while (! validate_z_ct0 (sig.z, ct0, sig.c, sk.t0));

            // Calculate (q/p)*xi2 = (q/p) * Round (c*s2) - c*s2, nu = Round ((p/q) * (xi1 - xi2))
            calculate_xi2_nu (nu, xi2, rcs2, sig.c, s2, xi1);

            // Calculate r1 and r0
            calculate_r (r, w, rcs2, nu);
            decompose (r1, r0, r, PQS_gamma_bar2_bits + 1);

        } while (! validate_r (r0, r1, w1));

        hw = make_hints (sig.h, ct0, r);

    } while (hw > PQS_omega);
}

// sign message m using keypair (sk, pk) and return the signed message into sm
int sign_from_keybytes (unsigned char *sm, unsigned long long *smlen, 
                        unsigned char *m, unsigned long long mlen, 
                        unsigned char *sk, unsigned char *pk) 
{
    signat_t sig;
    sk_t sk_obj;
    vk_t vk;

    #ifdef MATRIX_A
    vk.A = matrix_A;
    #else
    fprintf (stderr, "MATRIX_A undefined. Public A must be defined when using byte arrays for keys.\n");
    return -1;
    #endif

    deserialize_secret (sk_obj, sk);
    deserialize_verifykey (vk, pk);

    sign (sk_obj, vk, m, (unsigned int) mlen, sig);

    *smlen = pack_signed_message (sm, sig, m, mlen);

    if (*smlen != (mlen + SIGNATUREBYTES))
        return -2;

    return 0;
}

// SUBFUNCTIONS

// To make known answer tests reproducible, this code generates y from seed = H(K||mu), mu = H(tr||m)
void gen_seeds (unsigned char *mu, unsigned char *seed, sk_t &sk, unsigned char *m, unsigned int mlen) {
    int len1 = mlen + CRHBYTES;           // length of input to first hash
    int len2 = KBYTES + CRHBYTES;         // length of input to second hash

    unsigned char buf1[len1];
    unsigned char buf2[len2];

    // Calculate mu = hash (tr || m) 
    memcpy (buf1, sk.tr, CRHBYTES);
    memcpy (buf1 + CRHBYTES, m, mlen);
    shake256 (mu, CRHBYTES, buf1, mlen + CRHBYTES);

    // Calculate seed = hash (K || mu)
    memcpy (buf2, sk.K, KBYTES);
    memcpy (buf2 + KBYTES, mu, CRHBYTES);
    shake256 (seed, CRHBYTES, buf2, KBYTES + CRHBYTES);
}

void calculate_s2 (polyveck &s2, polyveck t, polymatkl A, polyvecl s1) {
    polyveck v;
    int32_t i, j, value;

    mult_MatVecl (v, A, s1, PQS_q);
    for (i = 0; i < PQS_k; i++) {
        for (j = 0; j < PQS_n; j++) {
            value = t.polynomial[i].coeffs[j] << PQS_qdivp_bits;
            s2.polynomial[i].coeffs[j] = value - v.polynomial[i].coeffs[j];
        }
    }
}

void generate_y (polyvecl &y, unsigned char *seed) {
    genVecl (y, PQS_gamma1 - 1, BUFLEN_gamma, seed);
}

void calculate_w_xi1 (polyveck &w, polyveck &xi1, polymatkl A, polyvecl y) {
    polyveck v;
    int32_t i, j, ayij, wij;

    mult_MatVecl (v, A, y, PQS_q);

    for (i = 0; i < PQS_k; i++) {
        for (j = 0; j < PQS_n; j++) {
            ayij = v.polynomial[i].coeffs[j];
            wij = w.polynomial[i].coeffs[j] = round (ayij);
            xi1.polynomial[i].coeffs[j] = (wij << PQS_qdivp_bits) - ayij;
        }
    }
}

void calculate_w1 (polyveck &w1, polyveck w) {
    int numbits = PQS_p_bits - (PQS_gamma_bar2_bits + 1);
    high_bits (w1, w, numbits);
}

void calculate_z (polyvecl &z, polyvecl y, polyvecl s1, poly c) {
    uint32_t i;

    for (i=0; i<PQS_l; i++) {
        MULT_POLYPOLY (c, s1.polynomial[i], z.polynomial[i], PQS_n, PQS_q);
        add_PolyPoly (z.polynomial[i], z.polynomial[i], y.polynomial[i], PQS_q);
    }
}

int validate_z_ct0 (polyvecl z, polyveck ct0, poly c, polyveck t0) {
    int32_t i, j, bound1, bound2, coeff;

    bound1 = PQS_gamma1 - PQS_beta1;
    bound2 = PQS_gamma_bar2;

    for (i = 0; i < PQS_k; i++) {
        MULT_POLYPOLY (c, t0.polynomial[i], ct0.polynomial[i], PQS_n, PQS_q);

        for (j = 0; j < PQS_n; j++) {
            coeff = z.polynomial[i].coeffs[j];
            if (coeff >= bound1 || -coeff <= -bound1) {

                printf ("failed condition 1: z[%d][%d] = %d\n", i, j, z.polynomial[i].coeffs[j]);
                return 0;
            }

            coeff = ct0.polynomial[i].coeffs[j];
            if (coeff >= bound2 || -coeff <= -bound2) {
                printf ("failed condition 2: ct0[%d][%d] = %d\n", i, j, ct0.polynomial[i].coeffs[j]);
                return 0;
            }
        }
    }

    return 1;
}

void calculate_xi2_nu (polyveck &nu, polyveck &xi2, polyveck &rcs2, poly c, polyveck s2, polyveck xi1) {
    int32_t i, j, cs2ij, rcs2ij, xi2ij;
    polyveck cs2;

    for (i = 0; i < PQS_k; i++) {
        MULT_POLYPOLY (c, s2.polynomial[i], cs2.polynomial[i], PQS_n, PQS_q);

        for (j = 0; j < PQS_n; j++) {
            cs2ij = cs2.polynomial[i].coeffs[j];
            rcs2ij = rcs2.polynomial[i].coeffs[j] = round (cs2ij);
            xi2ij = xi2.polynomial[i].coeffs[j] = (rcs2ij << PQS_qdivp_bits) - cs2ij;
            nu.polynomial[i].coeffs[j] = round (xi1.polynomial[i].coeffs[j] - xi2ij); 
        }
    }
}

void calculate_r (polyveck &r, polyveck w, polyveck rcs2, polyveck nu) {
    int i;

    for (i = 0; i < PQS_k; i++) {
        subt_PolyPoly (r.polynomial[i], w.polynomial[i], rcs2.polynomial[i], PQS_p); 
        subt_PolyPoly (r.polynomial[i], r.polynomial[i], nu.polynomial[i], PQS_p);
    }
}

int validate_r (polyveck r0, polyveck r1, polyveck w1) {
    int32_t i, j, bound, coeff;

    bound = PQS_gamma_bar2 - PQS_beta2;

    for (i = 0; i < PQS_k; i++) {
        for (j = 0; j < PQS_n; j++) {
            coeff = r0.polynomial[i].coeffs[j];
            if (coeff >= bound) {
                printf ("Failed condition 3: r0[%d][%d] = %d\n", i, j, coeff);
                return 0;
            }

            if (r1.polynomial[i].coeffs[j] != w1.polynomial[i].coeffs[j]) {
                printf ("Failed condition 4: r1[%d][%d] = %d, w1[%d][%d] = %d\n", i, j, r1.polynomial[i].coeffs[j], i, j, w1.polynomial[i].coeffs[j]);
                return 0;
            }
        }
    }

    return 1;
}

int make_hints (polyveck &h, polyveck ct0, polyveck r) {
    // returns hamming weight of hint vector
    uint32_t i, j, r1ij, v1ij, numbits, hw = 0;
    polyveck v, r1, v1;

    numbits = PQS_p_bits - (PQS_gamma_bar2_bits + 1);

    for (i = 0; i < PQS_k; i++)
        add_PolyPoly (v.polynomial[i], r.polynomial[i], ct0.polynomial[i], PQS_q);

    high_bits (r1, r, numbits);
    high_bits (v1, v, numbits);

    for (i = 0; i < PQS_k; i++) {
        for (j = 0; j < PQS_n; j++) {
            r1ij = r1.polynomial[i].coeffs[j];
            v1ij = v1.polynomial[i].coeffs[j];
            h.polynomial[i].coeffs[j] = (r1ij != v1ij);
            if (r1ij != v1ij)
                hw++;
        }
    }

    return hw;
}
