#ifndef PVPSVMC_BENCH_H
#define PVPSVMC_BENCH_H

#include "SVM.h"

void test_SVM_Gen();

void bench_ModelEnc(int features);

void bench_ModelEnc_improved(int features);

void load_data_poly(ModelPara &modelPara, vec_ZZ_p &y, Vec<vec_ZZ_p> &Z, int features, int slots);

void load_data_rbf(ModelPara &modelPara, vec_ZZ_p &y, Vec<vec_ZZ_p> &Z, int features, int slots);

void bench_compute_basic_poly(int features);

void bench_compute_basic_rbf(int features);

void bench_compute_improved(int features);

void bench_communication();



#endif //PVPSVMC_BENCH_H
