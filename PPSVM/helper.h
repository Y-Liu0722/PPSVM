#pragma once
#include <sstream>
#include <iostream>
#include <utility>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <ctime>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/BasicThreadPool.h>
#include <gmp.h>
#include <helib/helib.h>
#include <helib/matmul.h>

extern "C"{
    #include <relic/relic.h>
}
#define PI 3.141592654


using namespace std;
using namespace NTL;



void Random_ZZ_pX(ZZ_pX &a, int N, int q_bit);
void SecretKey(ZZ_pX &sk, int N, int hsk);
void GaussRand(ZZ_pX &e, int N);
Vec<vec_ZZ_p> GenRandomMatrix(int m, int n, const ZZ& r);

void DataProcess(double &mean, double &stdev, double *Time, int cyctimes);
void LCM(ZZ &res, ZZ a, ZZ b);
ZZ PRF_ZZ(int prfkey, ZZ mmod);
void ZZ2bn(bn_t out, ZZ in);

// Read data from csv file.
void read_csv(std::vector<std::vector<std::string>> &data, std::string filename);
void data_preprocess(std::vector<std::vector<double>> &f_data, std::vector<std::vector<std::string>> data);
void f_data2ZZ_p(Vec<vec_ZZ_p> &data, const vector<std::vector<double>>& f_data, double scale);

void bigEndianToHexString(uint8_t *data, size_t length, char *hexString);
void hexStringToBigEndian(const char *hexString, uint8_t *data, size_t length);
