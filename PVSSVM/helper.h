#pragma once
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <ctime>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <gmp.h>

extern "C"{
#include <relic/relic.h>
}

using namespace NTL;

void eval_at(bn_t v, const bn_t *poly, const bn_t x, const bn_t mod, size_t n);

void make_random_shares(bn_t *shares, bn_t *x, const bn_t secret, const bn_t mod, int n, int t);

void PI(bn_t &v, bn_t *inputs, size_t len);

void div_mod(bn_t v,const bn_t num,const bn_t den,const bn_t mod);

void lagrange_interpolate(bn_t secret, bn_t v, const bn_t *x, const bn_t* y, const bn_t mod, size_t k);

void bn2ZZ(ZZ &a, bn_t b);

typedef struct {
    bn_t *alpha; // (y_j alpha_j) j \in SV
    bn_t ** sv; // support vectors
    bn_t b; // bias

    // poly kernel
    bn_t gamma;
    bn_t c;

    long SV; // num of support vectors
    long features; // num of features
    bn_t scale; // scale factor
} ModelPara;


// Read data from csv file.
void read_csv(std::vector<std::vector<std::string>> &data, std::string filename);
void data_preprocess(std::vector<std::vector<double>> &f_data, std::vector<std::vector<std::string>> data);
void f_data2bn(bn_t **data, std::vector<std::vector<double>> f_data, double scale);

void set_model_paras(ModelPara &modelPara, double gamma, double b, double c, std::string para_file);
void get_user_inputs(bn_t **&X, bn_t *&y, std::string in_file, size_t &size);

void bn_lag(bn_t *c, const bn_t *a, const bn_t b, size_t n);

void bn_evl(bn_t c, const bn_t *a, const bn_t x, const bn_t b, size_t n);