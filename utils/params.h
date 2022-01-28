#ifndef PARAMS_H
#define PARAMS_H

#include <cstdint>

// This uses the medium parameter set

#define PQS_n 256

#define PQS_p_bits 19
#define PQS_q_bits 23
#define PQS_p (1 << PQS_p_bits)
#define PQS_q (1 << PQS_q_bits)

#define PQS_qdivp_bits (PQS_q_bits - PQS_p_bits)

#define PQS_d 10

#define PQS_gamma1 524288
#define PQS_gamma1_bits 19

#define PQS_gamma_bar1 32768
#define PQS_gamma_bar1_bits 15

#define PQS_gamma2 262144
#define PQS_gamma2_bits 18

#define PQS_gamma_bar2 16384
#define PQS_gamma_bar2_bits 14

#define PQS_eta 8
#define PQS_s1_bits 5

#define PQS_k 4
#define PQS_l 3

#define PQS_omega 80

#define PQS_beta1 425
#define PQS_beta2 25

#define SEEDBYTES 32
#define CRHBYTES 48
#define KBYTES 32
#define POLW1_SIZE_PACKED ((PQS_n*4)/8)

#define MAXRANDOM 33554431

#define S1BYTES ((PQS_n*PQS_l*PQS_s1_bits) / 8)
#define T0BYTES ((PQS_n*PQS_k*PQS_d) / 8)
#define T1BYTES ((PQS_n*PQS_k*(PQS_p_bits-PQS_d)) / 8)

#define ZBYTES ((PQS_n*PQS_l*(PQS_gamma1_bits+1)) / 8)
#define HBYTES (PQS_k + PQS_omega)
#define CBYTES ((PQS_n*2)/8) 

#define SECRETKEYBYTES (KBYTES + CRHBYTES + S1BYTES + T0BYTES)
#define PUBLICKEYBYTES (T1BYTES + SEEDBYTES)
#define SIGNATUREBYTES (ZBYTES + HBYTES + CBYTES)

// Number of bytes needed for rejection sampling, precomputed from Hoeffding's inequality rounded up to nearest multiple of SHAKE_RATE
#define BUFLEN_A 840 
#define BUFLEN_gamma 840
#define BUFLEN_s 336 


//#define DEBUG
//#define REGEN_A

// A polynomial in R_q, represented by a vector of coefficients (x^0, x^1, ..., x^(PQS_n-1))
typedef struct {
  int32_t coeffs[PQS_n];
} poly;

// A vector of length PQS_l of polynomials
typedef struct {
  poly polynomial[PQS_l];
} polyvecl;

// A vector of length PQS_k of polynomials
typedef struct {
  poly polynomial[PQS_k];
} polyveck;

// Matrix of size PQS_k*PQS_l of polynomials
typedef struct {
  poly elem[PQS_k][PQS_l];
} polymatkl;

// Secret key
typedef struct {
  unsigned char K[KBYTES];
  unsigned char tr[CRHBYTES];
  polyvecl s1;
  polyveck t0;
} sk_t;

// Public key
typedef struct {
  polymatkl A;
  polyveck t1;
} vk_t;

// Signature
typedef struct {
  polyvecl z;
  polyveck h;
  poly c;
} signat_t;

#endif
