#include "BKS.h"


void BKS_Gen(BKS_EK &ek1, BKS_EK &ek2, PKE_PK &pk, PKE_Para &pkePara) {
    PKE_SK sk;
    PKE_Gen(pkePara, pk, sk);

    ek1.prf.SetLength(2);
    ek1.s_b.SetLength(2);
    ek2.prf.SetLength(2);
    ek2.s_b.SetLength(2);

    Random_ZZ_pX(ek1.prf[0], pkePara.N, pkePara.q_bit);
    Random_ZZ_pX(ek1.prf[1], pkePara.N, pkePara.q_bit);
    ek2.prf[0] = ek1.prf[0];
    ek2.prf[1] = ek1.prf[1];

    Random_ZZ_pX(ek1.s_b[0], pkePara.N, pkePara.q_bit);
    Random_ZZ_pX(ek1.s_b[1], pkePara.N, pkePara.q_bit);

    ek2.s_b[0] = sk[0] - ek1.s_b[0];
    ek2.s_b[1] = sk[1] - ek1.s_b[1];
}

void BKS_Enc(Ciphertext &C, const PKE_Para &pkePara, const PKE_PK &pk, const ZZ_pX &x) {
    PKE_OKDM(C, pkePara, pk, x);
}

void BKS_Load(MemoryV &t_x, int b, const BKS_EK &ek, const PKE_Para &pkePara, const Ciphertext &C_x) {
    t_x.SetLength(2);
    PKE_DDec(t_x, pkePara, ek.s_b, C_x);
    if (b == 1){
        t_x[0] = t_x[0] - ek.prf[0];
        t_x[1] = t_x[1] - ek.prf[1];
    } else {
        t_x[0] = t_x[0] + ek.prf[0];
        t_x[1] = t_x[1] + ek.prf[1];
    }
}

void BKS_ADD1(MemoryV &t_z, int b, const BKS_EK &ek, const MemoryV &t_x, const MemoryV &t_y) {
    t_z.SetLength(2);
    if (b == 1){
        t_z[0] = t_x[0] + t_y[0] - ek.prf[0];
        t_z[1] = t_x[1] + t_y[1] - ek.prf[1];
    } else {
        t_z[0] = t_x[0] + t_y[0] + ek.prf[0];
        t_z[1] = t_x[1] + t_y[1] + ek.prf[1];
    }
}



void BKS_ADD2(Ciphertext &C_z, int b, const BKS_EK &ek, const Ciphertext &C_x, const Ciphertext &C_y) {
    C_z.c_m.SetLength(2);
    C_z.c_ms.SetLength(2);
    C_z.c_m[0] = C_x.c_m[0] + C_y.c_m[0];
    C_z.c_m[1] = C_x.c_m[1] + C_y.c_m[1];
    C_z.c_ms[0] = C_x.c_ms[0] + C_y.c_ms[0];
    C_z.c_ms[1] = C_x.c_ms[1] + C_y.c_ms[1];
}

void BKS_Mult(MemoryV &t_xy, int b, const BKS_EK &ek, const PKE_Para& pkePara, const MemoryV &t_x, const Ciphertext &C_y) {
    PKE_DDec(t_xy, pkePara, t_x, C_y);
    if (b == 1){
        t_xy[0] = t_xy[0] - ek.prf[0];
        t_xy[1] = t_xy[1] - ek.prf[1];
    } else {
        t_xy[0] = t_xy[0] + ek.prf[0];
        t_xy[1] = t_xy[1] + ek.prf[1];
    }
}

ZZ_pX BKS_Output(int b, const BKS_EK &ek, const MemoryV &t, const PKE_Para &pkePara) {
    ZZ_pX v;
    ZZ_p t_r;
    for (int i = 0; i < pkePara.N; i++){
        GetCoeff(t_r, t[0], i);
        t_r = t_r * 1;
        SetCoeff(v, i, t_r);
    }
    return v;
}

