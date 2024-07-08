#include "helper.h"

void Random_ZZ_pX(ZZ_pX &a, int N, int q_bit) {
    ZZ_p coeff;
    for (int i = 0; i < N; i++) {
        conv(coeff, RandomBits_ZZ(q_bit));
        SetCoeff(a, i, coeff);
    }
}

void SecretKey(ZZ_pX &sk, int N, int hsk) {
    int interval = 0;
    interval = N / hsk;
    int index = rand() % interval;
    for (int i = 0; i < hsk; i++) {
        SetCoeff(sk, index, 1);
        index = index + rand() % interval;
    }

}

void GaussRand(ZZ_pX &e, int N) {
    double res_standard;
    int deviation = 8;
    int res;
    for (int i = 0; i < N; i++) {
        res_standard = sqrt(-2.0 * log(rand() / (RAND_MAX + 1.0))) * sin(2.0 * PI * rand() / (RAND_MAX + 1.0));
        res = res_standard * deviation;
        SetCoeff(e, i, res);
    }
}

Vec<vec_ZZ_p> GenRandomMatrix(int m, int n, const ZZ& r){
    ZZ_p::init(r);
    Vec<vec_ZZ_p> M;
    M.SetLength(m);
    for (int i = 0; i < m; i++){
        M[i].SetLength(n);
        for (int j = 0; j < n; j++){
            M[i][j] = random_ZZ_p();
        }
    }

    return M;
}

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

void f_data2ZZ_p(Vec<vec_ZZ_p> &data, const vector<std::vector<double>>& f_data, double scale){
    for (const auto& r: f_data){
        vec_ZZ_p rows;
        for (const auto& cell: r){
//            std::cout << cell << " ";
            auto tmp = cell * scale;
            auto int_value = static_cast<long long>(std::round(tmp));
            rows.append(ZZ_p(int_value));
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
