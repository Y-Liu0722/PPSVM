#pragma once

#include "BKS.h"

typedef struct {
    vec_ZZ_p alpha; // (y_j alpha_j) j \in SV
    Vec<vec_ZZ_p> X; // support vectors
    ZZ_p b; // bias

    // kernel
    ZZ_p gamma;
    ZZ_p c;
    ZZ p;

    vec_ZZ_pX hat_X; // encoded support vectors
    ZZ_pX hat_alpha; // encoded alpha
    ZZ_pX hat_b; // encoded bias
    ZZ_pX hat_gamma; // encoded gamma
    ZZ_pX hat_c; // encoded c

    int m; // num of support vectors
    int n; // num of features
    int n_;
    ZZ_p scale; // scale factor
}ModelPara;

typedef ZZ PVK;

typedef struct {
    BKS_EK bksEk;
    Ciphertext C_delta;
}EK;

typedef struct {
    PKE_Para pkePara;
    ZZ g, r, r_;
    int m; // num of support vectors
    int n; // num of features
}PubPara;

vec_ZZ_p HomMVMult(int b, const EK &ek, const PKE_Para &pkePara, const Vec<MemoryV> &t_M, const Ciphertext &C_v, int m, int n);

void SVM_Gen(PubPara &para, PKE_PK &pk, EK &ek1, EK &ek2, PVK &pvk);

void SVM_ModelEnc_basic(Vec<Ciphertext> &C, Ciphertext& C_d, Ciphertext &C_b,
                        PubPara &para, const PKE_PK &pk, const ModelPara &modelPara);

void SVM_ModelEnc_improved(Vec<Ciphertext> &C, Ciphertext &C_d, Ciphertext &C_b, Ciphertext &C_g, Ciphertext &C_c,
                            PubPara &para, const PKE_PK &pk, const ModelPara &modelPara);

void Compute_basic_poly(ZZ_p &y_1, ZZ_p &y_2, ZZ &g_phi_1, ZZ &g_phi_2, const EK &ek1, const EK &ek2, const PubPara &para,
                        const Vec<Ciphertext> &C, const Ciphertext &C_d, const Ciphertext &C_b, const ZZ_p &gamma,
                        const ZZ_p &c, const ZZ &p, const PKE_PK &pk, const PVK &pvk, const vec_ZZ_p &z,
                        std::vector<double> &Time);

void Compute_basic_rbf_offline(Vec<MemoryV> &t_X_1, Vec<MemoryV> &t_X_2, Vec<MemoryV> &t_X_delta_1, Vec<MemoryV> &t_X_delta_2,
                                MemoryV &t_b_1, MemoryV &t_b_2, MemoryV &t_d_1, MemoryV &t_d_2,
                                vec_ZZ_p &H_1, vec_ZZ_p &H_2, vec_ZZ_p &Eta_1, vec_ZZ_p &Eta_2,
                                const EK &ek1, const EK &ek2, const PubPara &para,
                                const Vec<Ciphertext> &C, const Ciphertext& C_d, const Ciphertext &C_b, const ZZ_p& gamma, // server inputs
                                const PKE_PK& pk, const PVK& pvk, const vec_ZZ_p& z,
                                std::vector<double> &Time);

void Compute_basic_rbf(ZZ_p& y_1, ZZ_p& y_2, ZZ &g_phi_1, ZZ &g_phi_2, const EK &ek1, const EK &ek2, const PubPara &para,
                        const Vec<Ciphertext> &C, const Ciphertext& C_d, const Ciphertext &C_b, const ZZ_p& gamma, // server inputs
                        const PKE_PK& pk, const PVK& pvk, const vec_ZZ_p& z,
                        std::vector<double> &Time);

void Compute_improved(ZZ_p& y_1, ZZ_p& y_2, ZZ &g_phi_1, ZZ &g_phi_2, const EK &ek1, const EK &ek2, const PubPara &para,
                        const Vec<Ciphertext> &C, const Ciphertext& C_d, const Ciphertext &C_b,
                        const Ciphertext& C_gamma, const Ciphertext& C_c, const ZZ& p,// server inputs
                        const PKE_PK& pk, const vec_ZZ_p& z,
                        std::vector<double> &Time);