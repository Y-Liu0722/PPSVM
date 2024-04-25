//
// Created by ly on 24-3-28.
//

#ifndef PVSSVM_SSEDC_H
#define PVSSVM_SSEDC_H

#include "BGN.h"

void SSEDC_KeyGen(BGN_PK &pk, BGN_SK &sk);
void SSEDC_FGen(BGN_CT *CT_W, BGN_CT *CT_W2, bn_t *W, const BGN_PK &pk, size_t n);
void SSEDC_SGen(bn_t ** shares, bn_t *index, g1_t pvk, bn_t *X, BGN_PK &pk, int k, int t, size_t n);
void SSEDC_Compute(BGN_CT &v, BGN_CT &o, BGN_PK &pk, bn_t *shares, BGN_CT *CT_W, BGN_CT *CT_W2, int index, size_t n);
void SSEDC_PubVer(g1_t pvk, BGN_CT *CT_V, BGN_CT *CT_O, BGN_PK &pk, BGN_SK &sk, bn_t *index, int k, int t);


#endif //PVSSVM_SSEDC_H
