#ifndef PVPSVMC_BKS_H
#define PVPSVMC_BKS_H

#include "PKE.h"

typedef struct BKS_EK{
    vec_ZZ_pX prf;
    MemoryV s_b;
} BKS_EK;

void BKS_Gen(BKS_EK &ek1, BKS_EK &ek2, PKE_PK &pk, PKE_Para &pkePara);

void BKS_Enc(Ciphertext &C, const PKE_Para &pkePara, const PKE_PK &pk, const ZZ_pX &x);

void BKS_Load(MemoryV &t_x, int b, const BKS_EK &ek, const PKE_Para& pkePara, const Ciphertext &C_x);

void BKS_ADD1(MemoryV &t_z, int b, const BKS_EK &ek, const MemoryV &t_x, const MemoryV &t_y);

void BKS_ADD2(Ciphertext &C_z, int b, const BKS_EK &ek, const Ciphertext &C_x, const Ciphertext &C_y);

void BKS_Mult(MemoryV &t_xy, int b, const BKS_EK &ek, const PKE_Para& pkePara, const MemoryV &t_x, const Ciphertext &C_y);

ZZ_pX BKS_Output(int b, const BKS_EK &ek, const MemoryV &t, const PKE_Para &pkePara);
#endif //PVPSVMC_BKS_H
