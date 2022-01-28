#include <cstdio>

#include "keygen.h"

void keygen (sk_t &sk, vk_t &vk) {
    polyveck t;
  
    // Generate 32 byte secret K 
    randombytes (sk.K, KBYTES);

    // Generate secret s1
    generate_s1 (sk.s1);

    // Generate public matrix
    #ifdef REGEN_A
    generate_A (vk.A);
    #else
    vk.A = matrix_A;
    #endif

    // Calculate t and decompose into high bits (t1) and low bits (t0)
    calculate_t (t, vk.A, sk.s1);
    decompose (vk.t1, sk.t0, t, PQS_d);
    
    // Calculate tr
    calculate_tr (sk.tr, vk.t1);
}

// Generate keypair as byte arrays
int keygen_bytes (unsigned char *pk, unsigned char *sk) {
    sk_t sk_obj;
    vk_t vk_obj;

    keygen(sk_obj, vk_obj);

    serialize_secret (sk_obj, sk);
    serialize_verifykey (vk_obj, pk);

    return 0;
}

// SUBFUNCTIONS

void generate_s1 (polyvecl &s1) {
    unsigned char seed[SEEDBYTES];
    randombytes (seed, SEEDBYTES);
    genVecl(s1, PQS_eta, BUFLEN_s, seed);
}

void calculate_t (polyveck &t, polymatkl A, polyvecl s1) {
    // Multiply A and s1
    mult_MatVecl(t, A, s1, PQS_q);

    // Round each coefficient in each polynomial
    uint32_t i, j;
    for (i=0; i<PQS_k; i++) {
        for (j=0; j<PQS_n; j++) {
            t.polynomial[i].coeffs[j] = round (t.polynomial[i].coeffs[j]);
        }
    }
}

void calculate_tr (unsigned char *tr, polyveck t1) {
    unsigned char t1_bytes[T1BYTES];
    serialize_polyveck (t1, t1_bytes, PQS_p_bits - PQS_d);
    shake256 (tr, CRHBYTES, t1_bytes, T1BYTES);
}

void generate_A (polymatkl &A) {
    // generation of A accounts for about half the runtime of keygen and sign

    unsigned char seed[SEEDBYTES];

    // if DEBUG is defined then write the matrix to a file as it's generated
    #ifdef DEBUG
    FILE *mat_file;
    mat_file = fopen ("matrix_A.h", "w");

    fprintf (mat_file, "#ifndef MATRIX_A\n");
    fprintf (mat_file, "#define MATRIX_A\n\n");
    #endif

    poly p;
    uint32_t h, i, n;
    for (h=0; h<PQS_k; h++) {
        for (i=0; i<PQS_l; i++) {
            #ifdef DEBUG
            fprintf (mat_file, "static poly poly_%i_%i[PQS_n] = { ", h, i);
            #endif
            randombytes (seed, SEEDBYTES);
            rej_sample (p, PQS_q-1, BUFLEN_A, seed);
            for (n=0; n<PQS_n; n++) {
                if (p.coeffs[n] < 0) {
                    p.coeffs[n] = -p.coeffs[n];
                }

                if (h==0 && i == 0) {
                    p.coeffs[0] |= 1;
                    if (n > 0) {
                        p.coeffs[n] &= 0xFFFFFFFE;
                    }
                }

                #ifdef DEBUG
                if (n != PQS_n-1) {
                    fprintf (mat_file, "%7i, ", p.coeffs[n]);
                } else {
                    fprintf (mat_file, "%7i", p.coeffs[n]);
                }
                #endif
            }

            A.elem[h][i] = p;

            #ifdef DEBUG
            fprintf (mat_file, " };\n");
            #endif
        }
        #ifdef DEBUG
        fprintf (mat_file, "\n");
        #endif
    }

    #ifdef DEBUG
    fprintf (mat_file, "static polymatkl matrix_A = {");
    for (h=0; h<PQS_k; h++) {
        for (i=0; i<PQS_l; i++) {
            if ((h != PQS_k-1) | (i != PQS_l-1)) {
                fprintf (mat_file, "*poly_%i_%i, ", h, i);
            } else {
                fprintf (mat_file, "*poly_%i_%i", h, i);
            }
        }
    }

    fprintf (mat_file, "};\n\n");
    fprintf (mat_file, "#endif");
    fclose (mat_file);
    #endif
}
