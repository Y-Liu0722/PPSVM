#pragma once
#include "DJ.h"

typedef struct {
    int s;
    DJ_PK DJpk1;
    DJ_PK DJpk2;
    ZZ ct_;
    int K;

//    // verify para
//    ZZ g1_ord_ZZ;
//    ZZ g2_ord_ZZ;
//    bn_t g1_ord;
//    bn_t g2_ord;
//    ep_t g1_gen;
//    ep2_t g2_gen;
//    fp12_t gT_gen;
//    fp12_t const_gT;
    ZZ g; // generator of G
    ZZ r; // order of G
    ZZ r_; // module of G
} PVHSSPK;

typedef struct{
    ZZ phi;
    ZZ phi_;
    ep_t g1_gen;
    ep2_t g2_gen;
    ZZ C_A;
} PVHSSEK;

typedef struct {
    ZZ s; // share of phi * x
    ZZ s_; // share of phi * nu * xs
} REG;

typedef ZZ PVHSSPVK;

void PVHSS_Gen(PVHSSPK &pk, PVHSSEK &ek0, PVHSSEK &ek1, PVHSSPVK &pvk);

void PVHSS_Enc(ZZ &I, PVHSSPK pk, ZZ x);

void PVHSS_Load(REG &reg, int b, PVHSSPK pk, PVHSSEK ek, ZZ ct_x, int &prf_key);
void PVHSS_ADD(REG &reg, int b, PVHSSPK pk, PVHSSEK ek, REG reg_x, REG reg_y, int &prf_key);
void PVHSS_cMult(REG &reg, int b, PVHSSPK pk, PVHSSEK ek, REG reg_x, long long c, int &prf_key);
void PVHSS_Mult(REG &reg, int b, PVHSSPK pk, PVHSSEK ek, REG reg_x, ZZ ct_y, int &prf_key);
void PVHSS_Output(ZZ &v, ZZ &g_tau, int b, PVHSSPK pk, PVHSSEK ek, REG reg_x, int &prf_key);
bool PVHSS_Ver(ZZ &v_ZZ, ZZ v0, ZZ v1, ZZ g_tau0, ZZ g_tau1, PVHSSPK pk, PVHSSPVK pvk);