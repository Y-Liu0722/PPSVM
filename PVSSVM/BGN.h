
#pragma once
#include "helper.h"

typedef struct {
    bn_t q; // order
    g1_t g; // generator
    g2_t h; // generator
    g1_t gs; // g^s
    g1_t g1;
    g1_t g2;
}BGN_PK;

typedef struct {
    bn_t s;
    g1_t PI;
}BGN_SK;

typedef struct {
    g1_t c1;
    g1_t c2;
}BGN_CT;

void BGN_Gen(BGN_PK &pk, BGN_SK &sk);

void BGN_Enc(BGN_CT &ct, BGN_PK pk, bn_t m);

void BGN_Dec(bn_t &m, BGN_SK sk, BGN_CT ct);

void Sim_Dec(bn_t &res, bn_t target);

