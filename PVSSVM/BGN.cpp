#include "BGN.h"

void BGN_Gen(BGN_PK &pk, BGN_SK &sk) {
    core_init();
    ep_curve_init();
    ep_param_set(B12_P381);
    ep2_curve_init();
    ep2_curve_set_twist(EP_MTYPE);
    bn_null(pk.q);
    bn_new(pk.q);

    g1_get_ord(pk.q);
    g1_get_gen(pk.g);
    g2_get_gen(pk.h);

    bn_new(sk.s);
    bn_rand_mod(sk.s, pk.q);

    g1_mul_gen(pk.gs, sk.s);
    g1_rand(pk.g1);
    g1_rand(pk.g2);

    g1_t term1, term2;
    g1_mul(term1, pk.g1, sk.s);
    g1_neg(term2, pk.g2);
    g1_add(sk.PI, term1, term2);
}

void BGN_Enc(BGN_CT &ct, BGN_PK pk, bn_t m) {
    bn_t r;
    bn_new(r);
    bn_rand_mod(r, pk.q);

    g1_t term1, term2;
    g1_mul(term1, pk.g1, m);
    g1_mul(term2, pk.g, r);
    g1_add(ct.c1, term1, term2);

    g1_mul(term1, pk.g2, m);
    g1_mul(term2, pk.gs, r);
    g1_add(ct.c2, term1, term2);
}

void BGN_Dec(bn_t &m, BGN_SK sk, BGN_CT ct) {
    g1_t PI_C, tmp, t;
    g1_mul(PI_C, ct.c1, sk.s);
    g1_neg(tmp, ct.c2);
    g1_add(PI_C, PI_C, tmp);

    g1_copy(t, sk.PI);
    for (auto i = 0; i < 2147483647; i++) {
        if (g1_cmp(t, PI_C) == RLC_EQ) {
            bn_new(m);
            bn_set_dig(m, i + 1);
//            bn_print(m);
            break;
        }
        g1_add(t, t, sk.PI);
        g1_norm(t, t);
    }
}

void Sim_Dec(bn_t &res, bn_t target){
    bn_t bound;
    bn_new(bound);
    bn_set_2b(bound, 256);

    // binary search
    bn_t low, high, mid;
    bn_new(low);
    bn_new(high);
    bn_new(mid);
    bn_zero(low);
    bn_copy(high, bound);
//    bn_zero(res);
    while (bn_cmp(low, high) == RLC_LT) {
        bn_add(mid, low, high);
        bn_hlv(mid, mid);
//        std::cout << "mid: ";
//        bn_print(mid);
        if (bn_cmp(mid, target) == RLC_EQ) {
            bn_copy(res, mid);
            break;
        } else if (bn_cmp(mid, target) == RLC_LT) {
            bn_copy(low, mid);
        } else {
            bn_copy(high, mid);
        }
    }
}
