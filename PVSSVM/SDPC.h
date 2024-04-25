#include "BGN.h"

void SDPC_KeyGen(BGN_PK &pk, BGN_SK &sk);

void SDPC_FGen(BGN_CT *CT_W, bn_t *W, const BGN_PK &pk, size_t n);

void SDPC_SGen(bn_t ** shares, bn_t *index, g1_t pvk, bn_t *X, BGN_PK &pk, int k, int t, size_t n);

void SDPC_Compute(BGN_CT &v, BGN_CT &o, bn_t *shares, BGN_CT *CT_W, int index, size_t n);

void SDPC_PubVer(g1_t pvk, BGN_CT *CT_V, BGN_CT *CT_O, BGN_PK &pk, BGN_SK &sk, bn_t *index, int k, int t);