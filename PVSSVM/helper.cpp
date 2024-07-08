#include "helper.h"

// eval poly at x
// poly is coefficient of x^i
void eval_at(bn_t v, const bn_t *poly, const bn_t x, const bn_t mod, size_t n){
    bn_zero(v);
    // bn_print(v);
    for (int i = n - 1; i >= 0; i--){
        bn_mul(v, v, x);
        bn_mod(v, v, mod);
        bn_add(v, v, poly[i]);
        bn_mod(v, v, mod);
    }
}

// make random shares
// n server t threshold
void make_random_shares(bn_t *shares, bn_t *x, const bn_t secret, const bn_t mod, int n, int t){
    bn_t *poly = RLC_ALLOCA(bn_t, t+1);
    bn_copy(poly[0], secret);
//    std::cout << "poly 0:";
//    bn_print(poly[0]);
//    bn_print(mod);
    bn_t mod_t;
    bn_copy(mod_t, mod);
//    bn_print(mod_t);
    for (auto i = 1; i <= t; i++){
        bn_rand_mod(poly[i], mod_t);
//        std::cout << "poly " << i << ":";
//        bn_print(poly[i]);
    }

    bn_free(mod_t);
    for (auto i = 0; i < n; i++){
        bn_set_dig(x[i], i+1);
        bn_evl(shares[i], poly, x[i], mod, t+1);
//        eval_at(shares[i], poly, x[i], mod, t+1);
    }
//    delete[] poly;
}

// product of inputs
void PI(bn_t &v, bn_t *inputs, size_t len){
    bn_copy(v, inputs[0]);
    for (auto i = 1; i < len; i++){
        bn_mul(v, v, inputs[i]);
    }
}

// Compute num / den modulo prime p
void div_mod(bn_t v,const bn_t num,const bn_t den,const bn_t mod){
    bn_t c, d, e;
    bn_new(c);
    bn_new(d);
    bn_new(e);
    bn_gcd_ext_basic(c, d, e, den, mod);
    bn_mul(v, num, d);
    bn_free(c);
    bn_free(d);
    bn_free(e);
}

//// Lagrange interpolation
//void lagrange_interpolate(bn_t *poly, const bn_t *shares, const bn_t mod, size_t t){
//    bn_t *tmp = RLC_ALLOCA(bn_t, t + 1);
//    if (t == 0) {
//        bn_zero(poly[0]);
//        return;
//    }
//
//    for (int i = 0; i <= t; i++) {
//		bn_null(tmp[i]);
//		bn_new(tmp[i]);
//	}
//
//    for (int i = 0; i < t; i++) {
//        bn_zero(tmp[0]);
//        if (i == 0) {
//		    bn_set_dig(tmp[1], 1);
//            bn_sub(poly[0], mod, shares[i]);
//        } else {
//            for (int j = 0; j <= i; j++) {
//                bn_copy(tmp[j + 1], poly[j]);
//            }
//            for (int j = 0; j <= i; j++) {
//                bn_mul(poly[j], poly[j], shares[i]);
//                bn_mod(poly[j], poly[j], mod);
//                bn_sub(poly[j], tmp[j], poly[j]);
//                bn_mod(poly[j], poly[j], mod);
//            }
//        }
//        bn_copy(poly[i + 1], tmp[i + 1]);
//    }
//}

void lagrange_interpolate(bn_t secret, bn_t v, const bn_t *x, const bn_t* y, const bn_t mod, size_t k){
    bn_t * nums = RLC_ALLOCA(bn_t, k);
    bn_t * dens = RLC_ALLOCA(bn_t, k);
    for (int i = 0; i < k; i++){
        bn_null(nums[i]);
        bn_null(dens[i]);
        bn_new(nums[i]);
        bn_new(dens[i]);
    }
    for (int i = 0; i < k; i++){
        bn_t* others = RLC_ALLOCA(bn_t, k-1);
        bn_t* x_o = RLC_ALLOCA(bn_t, k-1);
        bn_t* cur_o = RLC_ALLOCA(bn_t, k-1);
        bn_t cur;
        bn_new(cur);
        for (int j = 0; j < k; j++){
            bn_null(others[j]);
            bn_new(others[j]);
            if (j < i){
                bn_copy(others[j], x[j]);
            } else if (j == i){
                bn_copy(cur, x[j]);
            } else {
                bn_copy(others[j-1], x[j]);
            }
        }
        for (int j =0; j < k - 1; j++){
            bn_null(x_o[j]);
            bn_new(x_o[j]);
            bn_null(cur_o[j]);
            bn_new(cur_o[j]);
            bn_sub(x_o[j], v, others[j]);
            bn_sub(cur_o[j], cur, others[j]);
            bn_mod(x_o[j], x_o[j], mod);
            bn_mod(cur_o[j], cur_o[j], mod);
        }
        PI(nums[i], x_o, k-1);
        PI(dens[i], cur_o, k-1);
    }
    bn_t den;
    bn_new(den);
    PI(den, dens, k);
    bn_t sum;
    bn_new(sum);
    bn_zero(sum);
    bn_t tmp;
    bn_new(tmp);
    for (int i = 0; i < k; i++){
        bn_mul(tmp, nums[i], den);
        bn_mul(tmp, tmp, y[i]);
        bn_mod(tmp, tmp, mod);
        div_mod(tmp, tmp, dens[i], mod);
        bn_add(sum, sum, tmp);
    }
    bn_mod(sum, sum, mod);
    div_mod(secret, sum, den, mod);
    bn_add(secret, secret, mod);
    bn_mod(secret, secret, mod);
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

void data_preprocess(std::vector<std::vector<double>> &f_data, std::vector<std::vector<std::string>> data) {
    for (const auto& r: data){
        std::vector<double> rows;
        rows.reserve(r.size());
        for (const auto& cell: r){
            rows.push_back(std::stod(cell));
        }
        f_data.push_back(rows);
    }
}


void f_data2bn(bn_t **data, std::vector<std::vector<double>> f_data, double scale){
    bool neg;
    for (auto i = 0; i < f_data.size(); i++){
        for (auto j = 0; j < f_data[0].size(); j++){
            auto tmp = f_data[i][j] * scale;
            auto int_value = static_cast<long long>(std::round(tmp));
            neg = false;
            if (int_value < 0) {
                int_value = - int_value;
                neg = true;
            }
            bn_read_str(data[i][j], std::to_string(int_value).c_str(), floor(log10(int_value)) + 1, 10);
            if (neg) {
                bn_neg(data[i][j], data[i][j]);
            }
        }
    }
}

void set_model_paras(ModelPara &modelPara, double gamma, double b, double c, std::string para_file) {
    std::vector<std::vector<std::string>> para_data;

    read_csv(para_data, para_file);
    std::vector<std::vector<double>> f_data;
    data_preprocess(f_data, para_data);

    modelPara.SV = f_data.size();
    modelPara.features = f_data[0].size() - 1;

    auto scale = pow(10, 10);

    auto int_scale = static_cast<long long>(std::round(scale));

    bn_read_str(modelPara.scale, std::to_string(int_scale).c_str(), floor(log10(int_scale))+1, 10);

//    bn_print(modelPara.scale);

    modelPara.alpha = new bn_t[modelPara.SV];
    for (int i = 0; i < modelPara.SV; i++){
        bn_new(modelPara.alpha[i]);
    }
    modelPara.sv = new bn_t*[modelPara.SV];
    for (int i = 0; i < modelPara.SV; i++){
        modelPara.sv[i] = new bn_t[modelPara.features];
        for (auto j = 0; j < modelPara.features; j++){
            bn_new(modelPara.sv[i][j]);
        }
    }

    auto bn_data = new bn_t*[f_data.size()];
    for (int i = 0; i < f_data.size(); i++){
        bn_data[i] = new bn_t[f_data[0].size()];
        for (auto j = 0; j < f_data[0].size(); j++){
            bn_new(bn_data[i][j]);
        }
    }

    f_data2bn(bn_data, f_data, scale);
    for(auto i = 0; i < modelPara.SV; i++){
        bn_copy(modelPara.alpha[i], bn_data[i][0]);
        for (auto j = 1; j < modelPara.features + 1; j++){
            bn_copy(modelPara.sv[i][j-1], bn_data[i][j]);
        }
    }
    auto int_gamma = static_cast<long long>(std::round(scale * gamma));
    auto int_b = static_cast<long long>(std::round(scale * b));
    auto int_c = static_cast<long long>(std::round(scale * c));

    bn_read_str(modelPara.gamma, std::to_string(int_gamma).c_str(), floor(log10(int_gamma))+1, 10);
    bn_read_str(modelPara.b, std::to_string(-int_b).c_str(), floor(log10(-int_b))+1, 10);
    bn_neg(modelPara.b, modelPara.b);
    bn_read_str(modelPara.c, std::to_string(int_c).c_str(), floor(log10(int_c))+1, 10);

//
//    for (auto i = 0; i < f_data.size(); i++){
//        RLC_FREE(bn_data[i]);
//    }
//    RLC_FREE(bn_data);
    delete[] bn_data;
}

void get_user_inputs(bn_t ** &X, bn_t* &y, std::string in_file, size_t &size) {
    std::vector<std::vector<std::string>> input_data;

    read_csv(input_data, in_file);
    std::vector<std::vector<double>> f_data;
    data_preprocess(f_data, input_data);

    auto scale = pow(10, 10);
    size = f_data.size();

    y = new bn_t[size];
    for (auto i = 0; i < size; i++){
        bn_new(y[i]);
    }
    X = new bn_t*[size];
    for (auto i = 0; i < size; i++){
        X[i] = new bn_t[f_data[0].size() - 1];
        for (auto j = 0; j < f_data[0].size() -1; j++){
            bn_new(X[i][j]);
        }
    }

    auto bn_data = new bn_t*[f_data.size()];
    for (int i = 0; i < f_data.size(); i++){
        bn_data[i] = new bn_t[f_data[0].size()];
        for (auto j = 0; j < f_data[0].size(); j++){
            bn_new(bn_data[i][j]);
        }
    }


    f_data2bn(bn_data, f_data, scale);
    for(auto i = 0; i < size; i++){
        bn_copy(y[i], bn_data[i][0]);
        for (auto j = 1; j < f_data[0].size() ; j++){
            bn_copy(X[i][j-1], bn_data[i][j]);
        }
    }

    for (auto i = 0; i < size; i++){
        delete []bn_data[i];
    }
    delete[] bn_data;
}

void bn2ZZ(ZZ &a, bn_t b){
    auto size = bn_size_str(b, 10);
    auto str = new char[size];
    bn_write_str(str, size, b, 10);
//    std::cout << str << std::endl;
    a = conv<ZZ>(str);
}


// @b module
void bn_lag(bn_t *c, const bn_t *a, const bn_t b, size_t n) {
    int i, j;
    bn_t *t = RLC_ALLOCA(bn_t, n + 1);

    if (n == 0) {
        bn_zero(c[0]);
        return;
    }

    if (t == NULL) {

    }
    for (i = 0; i <= n; i++) {
        bn_null(t[i]);
        bn_new(t[i]);
    }

    for (i = 0; i < n; i++) {
        bn_zero(t[0]);
        if (i == 0) {
            bn_set_dig(t[1], 1);
            bn_sub(c[0], b, a[i]);
        } else {
            for (j = 0; j <= i; j++) {
                bn_copy(t[j + 1], c[j]);
            }
            for (j = 0; j <= i; j++) {
                bn_mul(c[j], c[j], a[i]);
                bn_mod(c[j], c[j], b);
                bn_sub(c[j], t[j], c[j]);
                bn_mod(c[j], c[j], b);
            }
        }
        bn_copy(c[i + 1], t[i + 1]);
    }
}

void bn_evl(bn_t c, const bn_t *a, const bn_t x, const bn_t b, size_t n) {
    bn_zero(c);
    for (int j = n - 1; j >= 0; j--) {
        bn_mul(c, c, x);
        bn_mod(c, c, b);
        bn_add(c, c, a[j]);
        bn_mod(c, c, b);
    }
}
