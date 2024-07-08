#ifndef PVPSVMC_PKE_H
#define PVPSVMC_PKE_H

#include "helper.h"


typedef struct PKE_Para {
    int N,msg_bit,p_bit,q_bit;
    ZZ r,p,q,coeff,twice_p,twice_q,half_p;
    ZZ_pX xN;
    ZZ_pXModulus modulus, modulus_r;
    Vec<ZZ_pXModulus> modulus_f;
    ZZ_pContext r_context, q_context;
    vec_ZZ_pX factors; // length N
    vec_ZZ_pX e;
    ZZ_pX q_div_p;
    int hsk = 64;
//    int d;
//    int num_slots;
}PKE_Para;

typedef vec_ZZ_pX PKE_PK;
typedef vec_ZZ_pX PKE_SK;
typedef vec_ZZ_pX PKE_Ciphertext;
typedef struct Ciphertext {
    PKE_Ciphertext c_m;
    PKE_Ciphertext c_ms;
} Ciphertext;

typedef vec_ZZ_pX MemoryV;

//eval_poly=1, eval a poly
void SetPara(PKE_Para & pkePara);

void PKE_Gen(PKE_Para &pkePara, PKE_PK &pkePk, PKE_SK &pkeSk);

void PKE_Enc(PKE_Ciphertext &c, const PKE_Para& pkePara, const PKE_PK &pkePk, const ZZ_pX &m);

void PKE_OKDM(Ciphertext &C, const PKE_Para &pkePara, const PKE_PK &pkePk, const ZZ_pX &m);

void PKE_DDec(MemoryV &t_xy, const PKE_Para& pkePara, MemoryV t_x, Ciphertext C_y);

ZZ_pX Encode(const PKE_Para& pkePara, const vec_ZZ_p& x);

vec_ZZ_p Decode(const PKE_Para& pkePara, const ZZ_pX& hat_x);

#endif //PVPSVMC_PKE_H
