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

using namespace std;
using namespace NTL;

void DataProcess(double &mean, double &stdev, double *Time, int cyctimes);
void LCM(ZZ &res, ZZ a, ZZ b);
ZZ PRF_ZZ(int prfkey, ZZ mmod);
void ZZ2bn(bn_t out, ZZ in);

// Read data from csv file.
void read_csv(std::vector<std::vector<std::string>> &data, std::string filename);
void data_preprocess(std::vector<std::vector<double>> &f_data, std::vector<std::vector<std::string>> data);
void f_data2ZZ(Vec<Vec<ZZ>> &data, vector<std::vector<double>> f_data, double scale);

void bigEndianToHexString(uint8_t *data, size_t length, char *hexString);
void hexStringToBigEndian(const char *hexString, uint8_t *data, size_t length);

void FastPowerMod(ZZ &res, ZZ a, ZZ b, ZZ m);