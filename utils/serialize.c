#include "serialize.h"
#include <cstring>
#include <cstdlib>
#include <cstdio>

// Serialize and deserialize functions.
// There is no way to verify that arrays are large enough so it is up to the calling function to pass in sufficiently sized arrays.

void serialize_secret (sk_t &sk, unsigned char *sk_bytes) {
    memcpy (sk_bytes, sk.K, KBYTES);
    memcpy (sk_bytes + KBYTES, sk.tr, CRHBYTES);
    serialize_polyvecl (sk.s1, sk_bytes + KBYTES + CRHBYTES, PQS_s1_bits);
    serialize_polyveck (sk.t0, sk_bytes + KBYTES + CRHBYTES + S1BYTES, PQS_d);
}

void deserialize_secret (sk_t &sk, unsigned char *sk_bytes) {
    memcpy (sk.K, sk_bytes, KBYTES);
    memcpy (sk.tr, sk_bytes + KBYTES, CRHBYTES);
    deserialize_polyvecl (sk.s1, sk_bytes + KBYTES + CRHBYTES, PQS_s1_bits);
    deserialize_polyveck (sk.t0, sk_bytes + KBYTES + CRHBYTES + S1BYTES, PQS_d);
}

void serialize_verifykey (vk_t &vk, unsigned char *vk_bytes) {
    serialize_polyveck (vk.t1, vk_bytes, PQS_p_bits-PQS_d);
}

void deserialize_verifykey (vk_t &vk, unsigned char *vk_bytes) {
    deserialize_polyveck (vk.t1, vk_bytes, PQS_p_bits-PQS_d);
}

unsigned long long pack_signed_message (unsigned char *sm, signat_t &sig, unsigned char *m, unsigned long long mlen) {
    int i, j, num_nonzero, offset = mlen;

    // copy the message into sm
    memset (sm, 0, mlen + SIGNATUREBYTES);
    if (mlen != 0)
        memcpy (sm, m, (size_t) mlen);

    // insert sig.z into sm
    serialize_polyvecl (sig.z, sm + offset, PQS_gamma1_bits+1);
    offset += ZBYTES;

    // pack hints as k bytes: number of nonzero per polynomial + omega bytes: indices of nonzero coeffs)
    for (i = 0; i < PQS_k; i++) {
        num_nonzero = 0;

        for (j = 0; j < PQS_n; j++) {
            if (sig.h.polynomial[i].coeffs[j] != 0) {
                sm[offset+PQS_k+num_nonzero] = j;
                num_nonzero++;
            }
        }

        sm[offset+i] = num_nonzero;
    }

    // insert sig.c into sm
    offset += HBYTES;
    serialize_poly (sig.c, sm + offset, 2);

    return mlen + SIGNATUREBYTES;
}

unsigned long long unpack_signed_message (unsigned char *sm, signat_t &sig, unsigned char *m, unsigned long long smlen) {
    unsigned long long mlen = smlen - SIGNATUREBYTES;
    int i, j, index, num_nonzero, start = 0, offset = mlen;

    // unpack message
    memset (m, 0, (size_t) mlen);
    memcpy (m, sm, (size_t) mlen);

    // unpack sig.z
    deserialize_polyvecl (sig.z, sm + offset, PQS_gamma1_bits+1);
    offset += ZBYTES;

    // unpack hints
    for (i = 0; i < PQS_k; i++) {
        num_nonzero = sm[offset+i];

        for (j = start; j < num_nonzero; j++) {
            index = sm[offset+PQS_k+j];
            sig.h.polynomial[i].coeffs[index] = 1;
        }

        start = num_nonzero;
    }

    // unpack sig.c and convert back to range [-1,1]
    offset += HBYTES;
    deserialize_poly (sig.c, sm + offset, 2, true);

    return mlen;
}

int sig_compare (signat_t sig1, signat_t sig2) {
    int i;


    for (i = 0; i < PQS_l; i++)
        if (poly_compare(sig1.z.polynomial[i], sig2.z.polynomial[i], PQS_n) != 0)
            return 1;

    for (i = 0; i < PQS_k; i++)
        if (poly_compare(sig1.h.polynomial[i], sig2.h.polynomial[i], PQS_n) != 0)
            return 1;

    if (poly_compare(sig1.c, sig2.c, PQS_n) != 0)
        return 1;

    return 0;
}

int poly_compare (poly p1, poly p2, int length) {
    for (int i = 0; i < length; i++)
        if (p1.coeffs[i] != p2.coeffs[i]) {
            printf ("index %d: %d != %d\n", i, p1.coeffs[i], p2.coeffs[i]);
            return 1;
        }

    return 0;
}
