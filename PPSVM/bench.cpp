#include "bench.h"


void test_SVM_Gen(){
    PubPara para;
    PKE_PK pk;
    EK ek1, ek2;
    PVK pvk;

    para.pkePara.msg_bit = 32;

    // compute average time for running SVM_Gen 100 times
    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < 100; i++){
        SVM_Gen(para, pk, ek1, ek2, pvk);
    }

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Average time for running SVM_Gen 100 times: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 100 << "ms" << std::endl;
}


void load_data_rbf(ModelPara &modelPara, vec_ZZ_p &y, Vec<vec_ZZ_p> &Z, int features, int slots) {
    std::vector<std::vector<std::string>> data;
    read_csv(data, "../data/SMS_para_rbf_" + std::to_string(features) + ".csv");
    std::vector<std::vector<double>> f_data;
    data_preprocess(f_data, data);

    auto scale = pow(2, 13);
    Vec<vec_ZZ_p> ZZ_p_data;
    f_data2ZZ_p(ZZ_p_data, f_data, scale);
    modelPara.m = ZZ_p_data.length();
    modelPara.n = features;
    modelPara.n_ = ceil(modelPara.m / floor(slots / modelPara.n));
    vec_ZZ_p pad;
    pad.SetLength(features + 1, ZZ_p::zero());
    if (ZZ_p_data.length() < modelPara.n_ * int(floor(slots / modelPara.n))){
        // padding zero
        ZZ_p_data.SetLength(modelPara.n_ * int(floor(slots / modelPara.n)), pad);
    }

    for (auto row: ZZ_p_data){
        modelPara.alpha.append(row[0]);
        vec_ZZ_p row_t;
        for (int i = 1; i < row.length(); i++){
            row_t.append(row[i]);
        }
        modelPara.X.append(row_t);
    }

    switch (features) {
        case 64:
            // 1.2833947614433878
            // {'C': 1, 'coef0': 0, 'gamma': 'auto'}
            // bias: [-0.84353458]
            modelPara.gamma = ZZ_p(128);
            modelPara.b = ZZ_p(-6910);
            break;
        case 128:
            // 0.9782874329694846
            // {'C': 10, 'coef0': 0, 'gamma': 'scale'}
            // bias: [-0.69864086]
            modelPara.gamma = ZZ_p(64);
            modelPara.b = ZZ_p(-5723);
            break;
        case 256:
            // 0.6893312396389075
            // {'C': 10, 'coef0': 0, 'gamma': 'auto'}
            // bias: [-0.65250761]
            modelPara.gamma = ZZ_p(32);
            modelPara.b = ZZ_p(-5345);
            break;
        case 512:
            // 0.4690808856392363
            // {'C': 10, 'coef0': 0, 'gamma': 'auto'}
            // bias: [-0.52921213]
            modelPara.gamma = ZZ_p(16);
            modelPara.b = ZZ_p(-4335);
            break;
        case 1024:
            // 0.3038905750556731
            // {'C': 10, 'coef0': 0, 'gamma': 'auto'}
            // 0.9865470852017937
            // bias: [-0.15969503]
            modelPara.gamma = ZZ_p(8);
            modelPara.b = ZZ_p(-1308);
            break;
        default:
            std::cerr << "No match features !" << std::endl;
            return;
    }

    read_csv(data, "../data/SMS_test_data_rbf_" + std::to_string(features) + ".csv");
    data_preprocess(f_data, data);
    f_data2ZZ_p(ZZ_p_data, f_data, scale);

    for (auto row: ZZ_p_data){
        y.append(row[0]);
        vec_ZZ_p row_t;
        for (int i = 1; i < row.length(); i++){
            row_t.append(row[i]);
        }
        Z.append(row_t);
    }
}

void load_data_poly(ModelPara &modelPara, vec_ZZ_p &y, Vec<vec_ZZ_p> &Z, int features, int slots) {
//    ZZ r;
//    RandomPrime(r, 32);
//    while (r % (2 * 8192) != 1){
//        RandomPrime(r, 32);
//    }
//    ZZ_p::init(r);
    std::vector<std::vector<std::string>> data;
    read_csv(data, "../data/SMS_para_" + std::to_string(features) + ".csv");
    std::vector<std::vector<double>> f_data;
    data_preprocess(f_data, data);

    auto scale = pow(2, 13);
    Vec<vec_ZZ_p> ZZ_p_data;
    f_data2ZZ_p(ZZ_p_data, f_data, scale);
    modelPara.m = ZZ_p_data.length();
    modelPara.n = features;
    modelPara.n_ = ceil(modelPara.m / floor(slots / modelPara.n));
    vec_ZZ_p pad;
    pad.SetLength(features + 1, ZZ_p::zero());
    if (ZZ_p_data.length() < modelPara.n_ * int(floor(slots / modelPara.n))){
        // padding zero
        ZZ_p_data.SetLength(modelPara.n_ * int(floor(slots / modelPara.n)), pad);
    }

    for (auto row: ZZ_p_data){
        modelPara.alpha.append(row[0]);
        vec_ZZ_p row_t;
        for (int i = 1; i < row.length(); i++){
            row_t.append(row[i]);
        }
        modelPara.X.append(row_t);
    }

    switch (features) {
        case 64:
            // 1.2833947614433878
            // {'C': 0.1, 'coef0': 2, 'degree': 3, 'gamma': 'auto'}
            // bias: [-1.69750825]
            modelPara.gamma = ZZ_p(128);
            modelPara.c = ZZ_p(134217728);
            modelPara.p = ZZ(3);
            modelPara.b = ZZ_p(-13906);
            break;
        case 128:
            // 0.9782874329694846
            // {'C': 0.1, 'coef0': 2, 'degree': 3, 'gamma': 'scale'}
            // bias: [-1.52946633]
            modelPara.gamma = ZZ_p(64);
            modelPara.c = ZZ_p(134217728);
            modelPara.p = ZZ(3);
            modelPara.b = ZZ_p(-12529);
            break;
        case 256:
            // 0.6893312396389075
            // {'C': 1, 'coef0': 1, 'degree': 3, 'gamma': 'auto'}
            // bias: [-1.41767505]
            modelPara.gamma = ZZ_p(32);
            modelPara.c = ZZ_p(67108864);
            modelPara.p = ZZ(3);
            modelPara.b = ZZ_p(-11613);
            break;
        case 512:
            // 0.4690808856392363
            // {'C': 1, 'coef0': 1, 'degree': 3, 'gamma': 'auto'}
            // bias: [-1.35246086]
            modelPara.gamma = ZZ_p(16);
            modelPara.c = ZZ_p(67108864);
            modelPara.p = ZZ(3);
            modelPara.b = ZZ_p(-11709);
            break;
        case 1024:
            // 0.3038905750556731
            // {'C': 1, 'coef0': 2, 'degree': 3, 'gamma': 'auto'}
            // bias: [-1.39582207]
            modelPara.gamma = ZZ_p(8);
            modelPara.c = ZZ_p(134217728);
            modelPara.p = ZZ(3);
            modelPara.b = ZZ_p(-11434);
            break;
        default:
            std::cerr << "No match features !" << std::endl;
            return;
    }

    read_csv(data, "../data/SMS_test_data_" + std::to_string(features) + ".csv");
    data_preprocess(f_data, data);
    f_data2ZZ_p(ZZ_p_data, f_data, scale);

    for (auto row: ZZ_p_data){
        y.append(row[0]);
        vec_ZZ_p row_t;
        for (int i = 1; i < row.length(); i++){
            row_t.append(row[i]);
        }
        Z.append(row_t);
    }
}

void bench_ModelEnc(int features){
    PubPara para;
    PKE_PK pk;
    EK ek1, ek2;
    PVK pvk;

    para.pkePara.msg_bit = 32;

    SVM_Gen(para, pk, ek1, ek2, pvk);

    ModelPara modelPara;
    vec_ZZ_p y;
    Vec<vec_ZZ_p> Z;

    {
        ZZ_pPush push(para.pkePara.r_context);

        load_data_rbf(modelPara, y, Z, features, para.pkePara.N);

        Vec<vec_ZZ_p> X_;
        auto enc_time = GetTime();
        X_.SetLength(modelPara.n_);
        modelPara.hat_X.SetLength(modelPara.n_);
        for (auto i = 0; i < modelPara.n_; i++) {
//        X_[i].SetLength(para.pkePara.N);
            for (auto j = 0; j < ceil(para.pkePara.N / modelPara.n); j++) {
                X_[i].append(modelPara.X[i * para.pkePara.N / modelPara.n + j]);
            }
            modelPara.hat_X[i] = Encode(para.pkePara, X_[i]);
        }
        modelPara.hat_alpha = Encode(para.pkePara, modelPara.alpha);
        vec_ZZ_p b_;
        b_.SetLength(para.pkePara.N, ZZ_p::zero());
        b_[0] = modelPara.b;
        modelPara.hat_b = Encode(para.pkePara, b_);
        enc_time = GetTime() - enc_time;
        std::cout << "Feature:" << features << std::endl;
        std::cout << "Enc time: " << enc_time * 1000 << std::endl;
    }
    Vec<Ciphertext> C;
    Ciphertext C_d, C_b;

    std::vector<double> Time;
    double time;
    for (auto i = 0; i < 100; i++){
        time = GetTime();
        SVM_ModelEnc_basic(C, C_d, C_b, para, pk, modelPara);
        time = GetTime() - time;
        Time.push_back(time);
    }

    double mean, stdev;
    // mean of time
    mean = std::accumulate(Time.begin(), Time.end(), 0.0) / Time.size();
    // stdev of time
    stdev = std::sqrt(std::inner_product(Time.begin(), Time.end(), Time.begin(), 0.0) / Time.size() - mean * mean);

    std::cout << "Mean of time: " << mean*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time: " << stdev << std::endl;
}

void bench_ModelEnc_improved(int features) {
    PubPara para;
    PKE_PK pk;
    EK ek1, ek2;
    PVK pvk;

    para.pkePara.msg_bit = 128;

    SVM_Gen(para, pk, ek1, ek2, pvk);

    ModelPara modelPara;
    vec_ZZ_p y;
    Vec<vec_ZZ_p> Z;
    double enc_time;
    {
        ZZ_pPush push(para.pkePara.r_context);

        load_data_poly(modelPara, y, Z, features, para.pkePara.N);

        Vec<vec_ZZ_p> X_;
//        modelPara.n_ = ceil((modelPara.m * modelPara.n) / para.pkePara.N);
        X_.SetLength(modelPara.n_);
        modelPara.hat_X.SetLength(modelPara.n_);
        enc_time = GetTime();
        for (auto i = 0; i < modelPara.n_; i++) {
//        X_[i].SetLength(para.pkePara.N);
            for (auto j = 0; j < ceil(para.pkePara.N / modelPara.n); j++) {
                X_[i].append(modelPara.X[i * para.pkePara.N / modelPara.n + j]);
            }
            modelPara.hat_X[i] = Encode(para.pkePara, X_[i]);
        }
        modelPara.hat_alpha = Encode(para.pkePara, modelPara.alpha);
        vec_ZZ_p b_, g_, c_;
        b_.SetLength(para.pkePara.N, ZZ_p::zero());
        g_.SetLength(para.pkePara.N, modelPara.gamma);
        c_.SetLength(para.pkePara.N, ZZ_p::zero());
        b_[0] = modelPara.b;
        c_[0] = modelPara.c;
        modelPara.hat_b = Encode(para.pkePara, b_);
        modelPara.hat_gamma = Encode(para.pkePara, g_);
        modelPara.hat_c = Encode(para.pkePara, c_);
        enc_time = GetTime() - enc_time;
    }
//    std::cout << "Hello" << std::endl;
    Vec<Ciphertext> C;
    Ciphertext C_d, C_b, C_g, C_c;

    std::vector<double> Time;
    double time;
    for (auto i = 0; i < 100; i++){
        time = GetTime();
        SVM_ModelEnc_improved(C, C_d, C_b, C_g, C_c, para, pk, modelPara);
        time = GetTime() - time;
        Time.push_back(time);
    }

    double mean, stdev;
    // mean of time
    mean = std::accumulate(Time.begin(), Time.end(), 0.0) / Time.size();
    // stdev of time
    stdev = std::sqrt(std::inner_product(Time.begin(), Time.end(), Time.begin(), 0.0) / Time.size() - mean * mean);

    std::cout << "Feature:" << features << std::endl;
    std::cout << "Mean of time: " << mean*1000 << "ms" << std::endl;
    std::cout << "Time add Enc: " << (mean + enc_time) * 1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time: " << stdev << std::endl;
}

void bench_compute_basic_poly(int features){
    PubPara para;
    PKE_PK pk;
    EK ek1, ek2;
    PVK pvk;

    para.pkePara.msg_bit = 32;

    SVM_Gen(para, pk, ek1, ek2, pvk);

    ModelPara modelPara;
    vec_ZZ_p y;
    Vec<vec_ZZ_p> Z;

    {
        ZZ_pPush push(para.pkePara.r_context);

        load_data_poly(modelPara, y, Z, features, para.pkePara.N);

        Vec<vec_ZZ_p> X_;
//        modelPara.n_ = ceil((modelPara.m * modelPara.n) / para.pkePara.N);
        X_.SetLength(modelPara.n_);
        modelPara.hat_X.SetLength(modelPara.n_);
        for (auto i = 0; i < modelPara.n_; i++) {
//        X_[i].SetLength(para.pkePara.N);
            for (auto j = 0; j < ceil(para.pkePara.N / modelPara.n); j++) {
                X_[i].append(modelPara.X[i * para.pkePara.N / modelPara.n + j]);
            }
            modelPara.hat_X[i] = Encode(para.pkePara, X_[i]);
        }
        modelPara.hat_alpha = Encode(para.pkePara, modelPara.alpha);
        vec_ZZ_p b_;
        b_.SetLength(para.pkePara.N, ZZ_p::zero());
        b_[0] = modelPara.b;
        modelPara.hat_b = Encode(para.pkePara, b_);
    }

    Vec<Ciphertext> C;
    Ciphertext C_d, C_b;

    SVM_ModelEnc_basic(C, C_d, C_b, para, pk, modelPara);

    std::vector<double> s1_offline;
    std::vector<double> s2_offline;
    std::vector<double> s1_online;
    std::vector<double> s2_online;
    std::vector<double> user_enc;
    std::vector<double> user_compute;

    ZZ_p y_1, y_2;
    ZZ g_phi_1, g_phi_2;
    for (auto i = 0; i < 100; i ++){
        // 0 time offline for server 1
        // 1 time offline for server 2
        // 2 time online for server 1
        // 3 time online for server 2
        // 4 time user enc
        // 5 time user compute
        std::vector<double> Time;
        Time.resize(6);
        Compute_basic_poly(y_1, y_2, g_phi_1, g_phi_2, ek1, ek2, para, C, C_d, C_b,
                            modelPara.gamma, modelPara.c, modelPara.p, pk, pvk, Z[i], Time);
        s1_offline.push_back(Time[0]);
        s2_offline.push_back(Time[1]);
        s1_online.push_back(Time[2]);
        s2_online.push_back(Time[3]);
        user_enc.push_back(Time[4]);
        user_compute.push_back(Time[5]);
    }
    std::cout << "Feature: " << features <<std::endl;
    // mean s1 offline
    double mean_s1_offline = std::accumulate(s1_offline.begin(), s1_offline.end(), 0.0) / s1_offline.size();
    // stdev s1 offline
    double stdev_s1_offline = std::sqrt(std::inner_product(s1_offline.begin(), s1_offline.end(), s1_offline.begin(), 0.0) / s1_offline.size() - mean_s1_offline * mean_s1_offline);
    std::cout << "Mean of time for s1 offline: " << mean_s1_offline*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s1 offline: " << stdev_s1_offline << std::endl;

    // mean s2 offline
    double mean_s2_offline = std::accumulate(s2_offline.begin(), s2_offline.end(), 0.0) / s2_offline.size();
    // stdev s2 offline
    double stdev_s2_offline = std::sqrt(std::inner_product(s2_offline.begin(), s2_offline.end(), s2_offline.begin(), 0.0) / s2_offline.size() - mean_s2_offline * mean_s2_offline);
    std::cout << "Mean of time for s2 offline: " << mean_s2_offline*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s2 offline: " << stdev_s2_offline << std::endl;

    // mean s1 online
    double mean_s1_online = std::accumulate(s1_online.begin(), s1_online.end(), 0.0) / s1_online.size();
    // stdev s1 online
    double stdev_s1_online = std::sqrt(std::inner_product(s1_online.begin(), s1_online.end(), s1_online.begin(), 0.0) / s1_online.size() - mean_s1_online * mean_s1_online);
    std::cout << "Mean of time for s1 online: " << mean_s1_online*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s1 online: " << stdev_s1_online << std::endl;

    // mean s2 online
    double mean_s2_online = std::accumulate(s2_online.begin(), s2_online.end(), 0.0) / s2_online.size();
    // stdev s2 online
    double stdev_s2_online = std::sqrt(std::inner_product(s2_online.begin(), s2_online.end(), s2_online.begin(), 0.0) / s2_online.size() - mean_s2_online * mean_s2_online);
    std::cout << "Mean of time for s2 online: " << mean_s2_online*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s2 online: " << stdev_s2_online << std::endl;

    // mean user enc
    double mean_user_enc = std::accumulate(user_enc.begin(), user_enc.end(), 0.0) / user_enc.size();
    // stdev user enc
    double stdev_user_enc = std::sqrt(std::inner_product(user_enc.begin(), user_enc.end(), user_enc.begin(), 0.0) / user_enc.size() - mean_user_enc * mean_user_enc);
    std::cout << "Mean of time for user enc: " << mean_user_enc*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for user enc: " << stdev_user_enc << std::endl;

    // mean user compute
    double mean_user_compute = std::accumulate(user_compute.begin(), user_compute.end(), 0.0) / user_compute.size();
    // stdev user compute
    double stdev_user_compute = std::sqrt(std::inner_product(user_compute.begin(), user_compute.end(), user_compute.begin(), 0.0) / user_compute.size() - mean_user_compute * mean_user_compute);
    std::cout << "Mean of time for user compute: " << mean_user_compute*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for user compute: " << stdev_user_compute << std::endl;


}

void bench_compute_basic_rbf(int features){
    PubPara para;
    PKE_PK pk;
    EK ek1, ek2;
    PVK pvk;

    para.pkePara.msg_bit = 32;

    SVM_Gen(para, pk, ek1, ek2, pvk);

    ModelPara modelPara;
    vec_ZZ_p y;
    Vec<vec_ZZ_p> Z;

    {
        ZZ_pPush push(para.pkePara.r_context);

        load_data_rbf(modelPara, y, Z, features, para.pkePara.N);

        Vec<vec_ZZ_p> X_;
//        modelPara.n_ = ceil((modelPara.m * modelPara.n) / para.pkePara.N);
        X_.SetLength(modelPara.n_);
        modelPara.hat_X.SetLength(modelPara.n_);
        for (auto i = 0; i < modelPara.n_; i++) {
//        X_[i].SetLength(para.pkePara.N);
            for (auto j = 0; j < ceil(para.pkePara.N / modelPara.n); j++) {
                X_[i].append(modelPara.X[i * para.pkePara.N / modelPara.n + j]);
            }
            modelPara.hat_X[i] = Encode(para.pkePara, X_[i]);
        }
        modelPara.hat_alpha = Encode(para.pkePara, modelPara.alpha);
        vec_ZZ_p b_;
        b_.SetLength(para.pkePara.N, ZZ_p::zero());
        b_[0] = modelPara.b;
        modelPara.hat_b = Encode(para.pkePara, b_);
    }

    Vec<Ciphertext> C;
    Ciphertext C_d, C_b;

    SVM_ModelEnc_basic(C, C_d, C_b, para, pk, modelPara);

    std::vector<double> s1_offline;
    std::vector<double> s2_offline;
    std::vector<double> s1_online;
    std::vector<double> s2_online;
    std::vector<double> user_enc;
    std::vector<double> user_compute;

    ZZ_p y_1, y_2;
    ZZ g_phi_1, g_phi_2;
    for (auto i = 0; i < 100; i ++){
        // 0 time offline for server 1
        // 1 time offline for server 2
        // 2 time online for server 1
        // 3 time online for server 2
        // 4 time user enc
        // 5 time user compute
        std::vector<double> Time;
        Time.resize(6);
        Compute_basic_rbf(y_1, y_2, g_phi_1, g_phi_2, ek1, ek2, para, C, C_d, C_b,
                            modelPara.gamma, pk, pvk, Z[i], Time);
        s1_offline.push_back(Time[0]);
        s2_offline.push_back(Time[1]);
        s1_online.push_back(Time[2]);
        s2_online.push_back(Time[3]);
        user_enc.push_back(Time[4]);
        user_compute.push_back(Time[5]);
    }
    std::cout << "Feature: " << features <<std::endl;
    // mean s1 offline
    double mean_s1_offline = std::accumulate(s1_offline.begin(), s1_offline.end(), 0.0) / s1_offline.size();
    // stdev s1 offline
    double stdev_s1_offline = std::sqrt(std::inner_product(s1_offline.begin(), s1_offline.end(), s1_offline.begin(), 0.0) / s1_offline.size() - mean_s1_offline * mean_s1_offline);
    std::cout << "Mean of time for s1 offline: " << mean_s1_offline*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s1 offline: " << stdev_s1_offline << std::endl;

    // mean s2 offline
    double mean_s2_offline = std::accumulate(s2_offline.begin(), s2_offline.end(), 0.0) / s2_offline.size();
    // stdev s2 offline
    double stdev_s2_offline = std::sqrt(std::inner_product(s2_offline.begin(), s2_offline.end(), s2_offline.begin(), 0.0) / s2_offline.size() - mean_s2_offline * mean_s2_offline);
    std::cout << "Mean of time for s2 offline: " << mean_s2_offline*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s2 offline: " << stdev_s2_offline << std::endl;

    // mean s1 online
    double mean_s1_online = std::accumulate(s1_online.begin(), s1_online.end(), 0.0) / s1_online.size();
    // stdev s1 online
    double stdev_s1_online = std::sqrt(std::inner_product(s1_online.begin(), s1_online.end(), s1_online.begin(), 0.0) / s1_online.size() - mean_s1_online * mean_s1_online);
    std::cout << "Mean of time for s1 online: " << mean_s1_online*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s1 online: " << stdev_s1_online << std::endl;

    // mean s2 online
    double mean_s2_online = std::accumulate(s2_online.begin(), s2_online.end(), 0.0) / s2_online.size();
    // stdev s2 online
    double stdev_s2_online = std::sqrt(std::inner_product(s2_online.begin(), s2_online.end(), s2_online.begin(), 0.0) / s2_online.size() - mean_s2_online * mean_s2_online);
    std::cout << "Mean of time for s2 online: " << mean_s2_online*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s2 online: " << stdev_s2_online << std::endl;

    // mean user enc
    double mean_user_enc = std::accumulate(user_enc.begin(), user_enc.end(), 0.0) / user_enc.size();
    // stdev user enc
    double stdev_user_enc = std::sqrt(std::inner_product(user_enc.begin(), user_enc.end(), user_enc.begin(), 0.0) / user_enc.size() - mean_user_enc * mean_user_enc);
    std::cout << "Mean of time for user enc: " << mean_user_enc*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for user enc: " << stdev_user_enc << std::endl;

    // mean user compute
    double mean_user_compute = std::accumulate(user_compute.begin(), user_compute.end(), 0.0) / user_compute.size();
    // stdev user compute
    double stdev_user_compute = std::sqrt(std::inner_product(user_compute.begin(), user_compute.end(), user_compute.begin(), 0.0) / user_compute.size() - mean_user_compute * mean_user_compute);
    std::cout << "Mean of time for user compute: " << mean_user_compute*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for user compute: " << stdev_user_compute << std::endl;


}

void bench_compute_improved(int features){
    PubPara para;
    PKE_PK pk;
    EK ek1, ek2;
    PVK pvk;

    para.pkePara.msg_bit = 128;

    SVM_Gen(para, pk, ek1, ek2, pvk);

    ModelPara modelPara;
    vec_ZZ_p y;
    Vec<vec_ZZ_p> Z;

    {
        ZZ_pPush push(para.pkePara.r_context);

        load_data_poly(modelPara, y, Z, features, para.pkePara.N);

        Vec<vec_ZZ_p> X_;
//        modelPara.n_ = ceil((modelPara.m * modelPara.n) / para.pkePara.N);
        X_.SetLength(modelPara.n_);
        modelPara.hat_X.SetLength(modelPara.n_);
        for (auto i = 0; i < modelPara.n_; i++) {
//        X_[i].SetLength(para.pkePara.N);
            for (auto j = 0; j < ceil(para.pkePara.N / modelPara.n); j++) {
                X_[i].append(modelPara.X[i * para.pkePara.N / modelPara.n + j]);
            }
            modelPara.hat_X[i] = Encode(para.pkePara, X_[i]);
        }
        modelPara.hat_alpha = Encode(para.pkePara, modelPara.alpha);
        vec_ZZ_p b_, g_, c_;
        b_.SetLength(para.pkePara.N, ZZ_p::zero());
        b_[0] = modelPara.b;
        modelPara.hat_b = Encode(para.pkePara, b_);
        g_.SetLength(para.pkePara.N, modelPara.gamma);
        c_.SetLength(para.pkePara.N, ZZ_p::zero());
        c_[0] = modelPara.c;
        modelPara.hat_gamma = Encode(para.pkePara, g_);
        modelPara.hat_c = Encode(para.pkePara, c_);
    }

    Vec<Ciphertext> C;
    Ciphertext C_d, C_b, C_g, C_c;

    SVM_ModelEnc_improved(C, C_d, C_b, C_g, C_c,  para, pk, modelPara);

    std::vector<double> s1_offline;
    std::vector<double> s2_offline;
    std::vector<double> s1_online;
    std::vector<double> s2_online;
    std::vector<double> user_enc;

    ZZ_p y_1, y_2;
    ZZ g_phi_1, g_phi_2;
    for (auto i = 0; i < 100; i ++){
        // 0 time offline for server 1
        // 1 time offline for server 2
        // 2 time online for server 1
        // 3 time online for server 2
        // 4 time user enc
        // 5 time user compute
        std::vector<double> Time;
        Time.resize(6);
        Compute_improved(y_1, y_2, g_phi_1, g_phi_2, ek1, ek2, para, C, C_d, C_b,
                            C_g, C_c, modelPara.p, pk, Z[i], Time);
        s1_offline.push_back(Time[0]);
        s2_offline.push_back(Time[1]);
        s1_online.push_back(Time[2]);
        s2_online.push_back(Time[3]);
        user_enc.push_back(Time[4]);
    }
    std::cout << "Feature: " << features <<std::endl;
    // mean s1 offline
    double mean_s1_offline = std::accumulate(s1_offline.begin(), s1_offline.end(), 0.0) / s1_offline.size();
    // stdev s1 offline
    double stdev_s1_offline = std::sqrt(std::inner_product(s1_offline.begin(), s1_offline.end(), s1_offline.begin(), 0.0) / s1_offline.size() - mean_s1_offline * mean_s1_offline);
    std::cout << "Mean of time for s1 offline: " << mean_s1_offline*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s1 offline: " << stdev_s1_offline << std::endl;

    // mean s2 offline
    double mean_s2_offline = std::accumulate(s2_offline.begin(), s2_offline.end(), 0.0) / s2_offline.size();
    // stdev s2 offline
    double stdev_s2_offline = std::sqrt(std::inner_product(s2_offline.begin(), s2_offline.end(), s2_offline.begin(), 0.0) / s2_offline.size() - mean_s2_offline * mean_s2_offline);
    std::cout << "Mean of time for s2 offline: " << mean_s2_offline*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s2 offline: " << stdev_s2_offline << std::endl;

    // mean s1 online
    double mean_s1_online = std::accumulate(s1_online.begin(), s1_online.end(), 0.0) / s1_online.size();
    // stdev s1 online
    double stdev_s1_online = std::sqrt(std::inner_product(s1_online.begin(), s1_online.end(), s1_online.begin(), 0.0) / s1_online.size() - mean_s1_online * mean_s1_online);
    std::cout << "Mean of time for s1 online: " << mean_s1_online*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s1 online: " << stdev_s1_online << std::endl;

    // mean s2 online
    double mean_s2_online = std::accumulate(s2_online.begin(), s2_online.end(), 0.0) / s2_online.size();
    // stdev s2 online
    double stdev_s2_online = std::sqrt(std::inner_product(s2_online.begin(), s2_online.end(), s2_online.begin(), 0.0) / s2_online.size() - mean_s2_online * mean_s2_online);
    std::cout << "Mean of time for s2 online: " << mean_s2_online*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s2 online: " << stdev_s2_online << std::endl;

    // mean user enc
    double mean_user_enc = std::accumulate(user_enc.begin(), user_enc.end(), 0.0) / user_enc.size();
    // stdev user enc
    double stdev_user_enc = std::sqrt(std::inner_product(user_enc.begin(), user_enc.end(), user_enc.begin(), 0.0) / user_enc.size() - mean_user_enc * mean_user_enc);
    std::cout << "Mean of time for user enc: " << mean_user_enc*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for user enc: " << stdev_user_enc << std::endl;
}

void bench_communication(){
    PubPara para;
    PKE_PK pk;
    EK ek1, ek2;
    PVK pvk;

    para.pkePara.msg_bit = 32;

    SVM_Gen(para, pk, ek1, ek2, pvk);

    ZZ_p f;
    random(f);
    std::cout << NumBits(rep(f)) <<std::endl;

}
