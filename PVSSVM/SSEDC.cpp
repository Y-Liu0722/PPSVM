#include "SSEDC.h"

void SSEDC_KeyGen(BGN_PK &pk, BGN_SK &sk) {
    BGN_Gen(pk, sk);
}

void SSEDC_FGen(BGN_CT *CT_W, BGN_CT *CT_W2, bn_t *W, const BGN_PK &pk, size_t n) {
    for (int i = 0; i < n; i ++){
        BGN_Enc(CT_W[i], pk, W[i]);
        bn_t tmp;
        bn_new(tmp);
        bn_mul(tmp, W[i], W[i]);
        BGN_Enc(CT_W2[i], pk, tmp);
    }
}

void SSEDC_SGen(bn_t ** shares, bn_t *index, g1_t pvk, bn_t *X, BGN_PK &pk, int k, int t, size_t n){
    bn_t alpha;
    bn_t mod;
    bn_new(mod);
    bn_set_dig(mod, 524287);
    bn_new(alpha);
    bn_rand_mod(alpha, mod);
    make_random_shares(shares[n], index, alpha, mod, k, t);
    for (int i = 0; i < n; i++){
        make_random_shares(shares[i], index, X[i], mod, k, t);
    }
    g1_mul_gen(pvk, alpha);
}

void SSEDC_Compute(BGN_CT &v, BGN_CT &o, BGN_PK &pk,bn_t *shares, BGN_CT *CT_W, BGN_CT *CT_W2, int index, size_t n) {
    g1_sub(v.c1, v.c1, v.c1);
    g1_sub(v.c2, v.c2, v.c2);
    g1_t t1, t2;
    g1_new(t1); g1_new(t2);
    bn_t tmp;
    bn_new(tmp);
    for (int i = 0; i < n; i++){
        g1_add(v.c1, v.c1, CT_W2[i].c1);
        g1_add(v.c2, v.c2, CT_W2[i].c2);

        bn_mul_dig(tmp, shares[i], 2);
        bn_neg(tmp, tmp);
        g1_mul(t1, CT_W[i].c1, tmp);
        g1_mul(t2, CT_W[i].c2, tmp);
        g1_add(v.c1, v.c1, t1);
        g1_add(v.c2, v.c2, t2);

        bn_mul(tmp, shares[i], shares[i]);
        g1_mul(t1, pk.g1, tmp);
        g1_mul(t2, pk.g2, tmp);

        g1_add(v.c1, v.c1, t1);
        g1_add(v.c2, v.c2, t2);
    }

    g1_norm(v.c1, v.c1);
    g1_norm(v.c2, v.c2);
    g1_free(t1); g1_free(t2);

    g1_mul(o.c1, v.c1, shares[n]);
    g1_mul(o.c2, v.c2, shares[n]);
    g1_norm(o.c1, o.c1);
    g1_norm(o.c2, o.c2);
}

void SSEDC_PubVer(g1_t pvk, BGN_CT *CT_V, BGN_CT *CT_O, BGN_PK &pk, BGN_SK &sk, bn_t *index, int k, int t){
    bn_t *V = RLC_ALLOCA(bn_t, k);
    bn_t *O = RLC_ALLOCA(bn_t, k);
    for (int i = 0; i < k; i++){
        bn_null(V[i]);
        bn_null(O[i]);
        bn_new(V[i]);
        bn_new(O[i]);
        BGN_Dec(V[i], sk, CT_V[i]);
        BGN_Dec(O[i], sk, CT_O[i]);
    }

    bn_t zero;
    bn_new(zero);
    bn_zero(zero);
    bn_t phi, psi;
    bn_new(phi);
    bn_new(psi);
    bn_t mod;
    bn_new(mod);
    bn_set_dig(mod, 524287);
    lagrange_interpolate(phi, zero, index, V, mod, t);
    lagrange_interpolate(psi, zero, index, O, mod, t);
//    bn_print(phi);

    g1_t left_side, right_side;
    g1_new(left_side); g1_new(right_side);
    g1_mul(left_side, pk.g, psi);
    g1_mul(right_side, pvk, phi);
    if (g1_cmp(left_side, right_side) == RLC_EQ){
        std::cout << "True" << std::endl;
    } else {
        std::cout << "Error" << std::endl;
    }

}