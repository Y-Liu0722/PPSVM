#pragma once

#include "PVHSS.h"

typedef struct {
    Vec<ZZ> alpha; // (y_j alpha_j) j \in SV
    Vec<Vec<ZZ>> sv; // support vectors
    ZZ b; // bias

    // poly kernel
    ZZ gamma;
    ZZ c;

    long SV; // num of support vectors
    long features; // num of features
    ZZ scale; // scale factor
}ModelPara;

void set_model_paras(ModelPara &modelPara, double gamma, double b, double c, std::string para_file);
void get_user_inputs(Vec<Vec<ZZ>> &X, Vec<ZZ> &y, std::string in_file);

void poly_kernel(REG &reg, int b, PVHSSPK pk, PVHSSEK ek,
                 Vec<ZZ> ct_x, ZZ ct_ccc, ZZ ct_3ccg,
                 ZZ ct_3cgg, ZZ ct_ggg, Vec<ZZ> ct_z,
                 int &prf_key, std::ofstream &bench_time);

void predict(REG &reg, int b, PVHSSPK pk, PVHSSEK ek, int &prf_key,
             Vec<Vec<ZZ>> ct_X, Vec<ZZ> ct_z, Vec<ZZ> ct_alpha,
             ZZ ct_ccc, ZZ ct_3ccg, ZZ ct_3cgg, ZZ ct_ggg, ZZ ct_b);

void eval_svm(std::string in_file, std::string para_file, double gamma, double b, double c);

void dot_prod(REG &reg, int b, PVHSSPK pk, PVHSSEK ek, int &prf_key, Vec<ZZ> ct_x, Vec<ZZ> ct_z, std::ofstream &dot_prod_time);

void squared_euclidean_distance(REG &reg, int b, PVHSSPK pk, PVHSSEK ek, int &prf_key, Vec<ZZ> ct_x, Vec<ZZ> ct_z, Vec<ZZ> ct_zz, std::ofstream &squared_euclidean_distance_time);


void test_input_dp_time(std::string in_file, std::string type, int features);

void test_input_sv_time(std::string para_file, std::string type, int features);

void test_batch_verify(int size);

void test_verify(int size);

void test_input_size();