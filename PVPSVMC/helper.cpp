#include "helper.h"

#include <utility>

void DataProcess(double &mean, double &stdev, double *Time, int cyctimes)
{
    double temp;
    double sum = 0;
    for (int i = 0; i < cyctimes; i++)
    {
        sum = sum + Time[i];
    }
    mean = sum / cyctimes;
    double temp_sum = 0;
    for (int i = 0; i < cyctimes; i++)
    {
        temp = mean - Time[i];
        temp = temp * temp;
        temp_sum = temp_sum + temp;
    }
    stdev = sqrt(temp_sum / cyctimes);
    stdev = stdev / mean;
}

ZZ PRF_ZZ(int prfkey, ZZ mmod) {
    ZZ res;
    SetSeed(ZZ(prfkey));
    RandomBnd(res, mmod);
    return res;
}

void LCM(ZZ &res, ZZ a, ZZ b)
{
    ZZ tmp;
    mul(tmp, a, b);
    div(res, tmp, GCD(a, b));
}



void read_csv(std::vector<std::vector<std::string>> &data, std::string filename) {
    std::ifstream file(filename);

    // check that the file open correctly
    if (!file.is_open()) {
        std::cerr << "Open file failed" << std::endl;
        return ;
    }

    // read csv file by rows
    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> rowData;
        std::istringstream iss(line);
        std::string cell;

        while (std::getline(iss, cell, ',')) {
            rowData.push_back(cell);
        }

        data.push_back(rowData);
    }

    file.close();
}

void data_preprocess(vector<std::vector<double>> &f_data, std::vector<std::vector<std::string>> data) {
    for (const auto& r: data){
        std::vector<double> rows;
        rows.reserve(r.size());
        for (const auto& cell: r){
            rows.push_back(std::stod(cell));
        }
        f_data.push_back(rows);
    }
}

void f_data2ZZ(Vec<Vec<ZZ>> &data, vector<std::vector<double>> f_data, double scale){
    for (const auto& r: f_data){
        Vec<ZZ> rows;
        for (const auto& cell: r){
//            std::cout << cell << " ";
            auto tmp = cell * scale;
            auto int_value = static_cast<long long>(std::round(tmp));
            rows.append(ZZ(int_value));
        }
//        std::cout << std::endl;
//        std::cout << rows << std::endl;
        data.append(rows);
    }
}

void bigEndianToHexString(uint8_t *data, size_t length, char *hexString) {
    for (size_t i = 0; i < length; i++) {
        sprintf(hexString + i * 2, "%02X", data[i]);
    }
}

void hexStringToBigEndian(const char *hexString, uint8_t *data, size_t length) {
    for (size_t i = 0; i < length; i++) {
        sscanf(hexString + i * 2, "%2hhX", &data[i]);
    }
}

void FastPowerMod(ZZ &res, ZZ a, ZZ b, ZZ m){
    res = ZZ(1);
    while (b > 0){
        if (b % 2 == 1){
            MulMod(res, res, a, m);
        }
        MulMod(a, a, a, m);
        b = b / 2;
    }
}
