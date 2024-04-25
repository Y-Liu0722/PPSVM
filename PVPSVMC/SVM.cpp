#include "SVM.h"

// read the model para from file and inputs
void set_model_paras(ModelPara &modelPara, double gamma, double b, double c, std::string para_file) {
    std::vector<std::vector<std::string>> para_data;

    read_csv(para_data, para_file);
    std::vector<std::vector<double>> f_data;
    data_preprocess(f_data, para_data);
    auto scale = pow(10, 15);
    power(modelPara.scale, 10, 15);

    Vec<Vec<ZZ>> ZZ_data;
    f_data2ZZ(ZZ_data, f_data, scale);

    for (auto r : ZZ_data){
        modelPara.alpha.append(r[0]);
        Vec<ZZ> tmp;
        for (auto i = 1; i < r.length(); i++){
            tmp.append(r[i]);
        }
        modelPara.sv.append(tmp);
    }

    modelPara.gamma = ZZ(static_cast<long long>(std::round(scale * gamma)));
    modelPara.b = ZZ(static_cast<long long>(std::round( scale * b)));
    mul(modelPara.b, modelPara.b, power(modelPara.scale, 9));
    modelPara.c = ZZ(static_cast<long long>(std::round(scale * c)));
    mul(modelPara.c, modelPara.c, power(modelPara.scale, 2));
    modelPara.SV = modelPara.sv.length();
    modelPara.features = modelPara.sv[0].length();
}


// compute kernel functions
// c^3 + 3c^2 \cdot \gamma \cdot \left(\sum_{i=1}^n x_{j,i} \cdot z_i\right) \\
// + 3c\cdot \gamma^2 \cdot \left(\sum_{i=1}^n \sum_{u=1}^n x_{j,i} \cdot z_i \cdot x_{j,u} \cdot z_u\right) \\
// + 3\gamma^3 \cdot \left(\sum_{i=1}^n \sum_{u=1}^n \sum_{k=1}^n  x_{j,i} \cdot z_j \cdot x_{j,u} \cdot z_u \cdot x_{j,k} \cdot z_k\right)
void poly_kernel(REG &reg, int b, PVHSSPK pk, PVHSSEK ek,
                 Vec<ZZ> ct_x, ZZ ct_ccc, ZZ ct_3ccg,
                 ZZ ct_3cgg, ZZ ct_ggg, Vec<ZZ> ct_z,
                 int &prf_key, std::ofstream &bench_time) {
    // Precompute all x_i and x_i * x_j and x_i * x_j * x_k
    Vec<REG> reg_x;
    long features = ct_x.length();
    reg_x.SetLength(features);
    for (auto i = 0; i < features; i++){
        if (ct_x[i] == 0){
            reg_x[i].s = 0;
            reg_x[i].s_ = 0;
            continue;
        }
        PVHSS_Load(reg_x[i], b, pk, ek, ct_x[i], prf_key);
    }

    Vec<Vec<REG>> reg_xx;
    reg_xx.SetLength(features);
    for (auto i = 0; i < features; i++){
        reg_xx[i].SetLength(features);
        for (auto j = 0; j < features; j++){
            if (ct_x[i] == 0 || ct_x[j] == 0){
                reg_xx[i][j].s = 0;
                reg_xx[i][j].s_ = 0;
                continue;
            }
            PVHSS_Mult(reg_xx[i][j], b, pk, ek, reg_x[i], ct_x[j], prf_key);
        }
    }

    Vec<Vec<REG>> reg_xz_xz;
    reg_xz_xz.SetLength(features);
    for (auto i = 0; i < features; i++){
        reg_xz_xz[i].SetLength(features);
    }

    // PreCompute c^3
    REG term1;
    PVHSS_Load(term1, b, pk, ek, ct_ccc, prf_key); // c^3

    // Start online computing
    std::cout << "Start online computing" << std::endl;
    double time = GetTime();

    // 3c^2 \cdot \gamma \cdot \left(\sum_{i=1}^n x_{j,i} \cdot z_i\right)
    REG term2;
    for (auto i = 0; i < features; i++){
        if ((reg_x[i].s == 0 && reg_x[i].s_ == 0) || ct_z[i] == 0){
            continue;
        }
        REG tmp;
        PVHSS_Mult(tmp, b, pk, ek, reg_x[i], ct_z[i], prf_key); // x_i * z_i
        PVHSS_ADD(term2, b, pk, ek, term2, tmp, prf_key); // term2 += x_i * z_i
    } // term2 = \sum x_i * z_i
    if (term2.s == 0 && term2.s_ == 0){
        PVHSS_ADD(reg, b, pk, ek, reg, term1, prf_key);
        time = GetTime() - time;
        std::cout << "All zero, early return" << std::endl;
        std::cout << "Online computing time: " << time * 1000 << " ms" << std::endl;
        bench_time << time * 1000 << std::endl;
        return;
    }

    PVHSS_Mult(term2, b, pk, ek, term2, ct_3ccg, prf_key); // 3c^2 \cdot \gamma \cdot * (\sum x_i * z_i)

    // 3c\cdot \gamma^2 \cdot \left(\sum_{i=1}^n \sum_{u=1}^n x_{j,i} \cdot z_i \cdot x_{j,u} \cdot z_u\right)
    REG term3;
    for (auto i = 0; i < features; i++){
        if ((reg_xx[i][i].s == 0 && reg_xx[i][i].s_ == 0) || ct_z[i] == 0){
            reg_xz_xz[i][i].s = 0;
            reg_xz_xz[i][i].s_ = 0;
            continue;
        }
        REG tmp;
        PVHSS_Mult(tmp, b, pk, ek, reg_xx[i][i], ct_z[i], prf_key);// x_i * x_i * z_i
        PVHSS_Mult(tmp, b, pk, ek, tmp, ct_z[i], prf_key);// x_i * z_i * x_i * z_i
        reg_xz_xz[i][i] = tmp; // reg_xz_xz[i][i] = x_i * z_i * x_i * z_i
        PVHSS_ADD(term3, b, pk, ek, term3, tmp, prf_key); // term3 += x_i * z_i * x_i * z_i
    } // \sum x_i * z_i * x_j * z_j when i=j
    for (auto i = 0; i <features; i++){
        for (auto j = i + 1; j < features; j++){
            if ((reg_xx[i][j].s == 0 && reg_xx[i][j].s_ == 0) || ct_z[i] == 0 || ct_z[j] == 0){
                reg_xz_xz[i][j].s = 0;
                reg_xz_xz[i][j].s_ = 0;
                continue;
            }
            REG tmp;
            PVHSS_Mult(tmp, b, pk, ek, reg_xx[i][j], ct_z[i], prf_key);// x_j * x_i * z_i
            PVHSS_Mult(tmp, b, pk, ek, tmp, ct_z[j], prf_key);// x_j * z_j * x_i * z_i
            reg_xz_xz[i][j] = tmp; // reg_xz_xz[i][i] = x_i * z_i * x_j * z_j
            // each  x_j * z_j * x_i * z_i appear twice
            PVHSS_cMult(tmp, b, pk, ek, tmp, 2, prf_key); // 2 * x_j * z_j * x_i * z_i
            PVHSS_ADD(term3, b, pk, ek, term3, tmp, prf_key); // term3 += 2 * x_j * z_j * x_i * z_i
        }
    }// term3 = \sum x_i * z_i * x_j * z_j
    PVHSS_Mult(term3, b, pk, ek, term3, ct_3cgg, prf_key); // 3c\cdot \gamma^2 * (\sum x_i * z_i * x_j * z_j)

    // 3\gamma^3 \cdot \left(\sum_{i=1}^n \sum_{u=1}^n \sum_{k=1}^n  x_{j,i} \cdot z_j \cdot x_{j,u} \cdot z_u \cdot x_{j,k} \cdot z_k\right)
    REG term4;
    for (auto i = 0; i < features; i++){
        if ((reg_xz_xz[i][i].s == 0 && reg_xz_xz[i][i].s_ == 0) || ct_x[i] == 0 || ct_z[i] == 0){
            continue;
        }
        REG tmp;
        PVHSS_Mult(tmp, b, pk, ek, reg_xz_xz[i][i], ct_z[i], prf_key);// x_i * x_i * z_i * z_i
        PVHSS_Mult(tmp, b, pk, ek, tmp, ct_x[i], prf_key);// x_i * z_i * x_i * z_i * x_i * z_i
        PVHSS_ADD(term4, b, pk, ek, term4, tmp, prf_key); // term4 += x_i * z_i * x_i * z_i * x_i * z_i
    } // \sum x_i * z_i * x_j * z_j * x_k * z_k when i=j=k
    for (auto i =0; i < features; i++){
        for (auto j = i + 1; j < features; j++){
            if ((reg_xz_xz[i][j].s == 0 && reg_xz_xz[i][j].s_ == 0) || ct_x[j] == 0 || ct_z[j] == 0){
                continue;
            }
            REG tmp;
            PVHSS_Mult(tmp, b, pk, ek, reg_xz_xz[i][i], ct_z[j], prf_key);// z_j * x_i * z_i * x_i * z_i
            PVHSS_Mult(tmp, b, pk, ek, tmp, ct_x[j], prf_key);// x_j * z_j * x_i * z_i * x_i * z_i
            // appear 3 times s.a. 112 121 211
            PVHSS_cMult(tmp, b, pk, ek, tmp, 3, prf_key); // 3 * x_j * z_j * x_i * z_i * x_i * z_i
            PVHSS_ADD(term4, b, pk, ek, term4, tmp, prf_key); // term4 += 3 * x_j * z_j * x_i * z_i * x_i * z_i

            PVHSS_Mult(tmp, b, pk, ek, reg_xz_xz[j][j], ct_z[i], prf_key);// z_i * x_j * z_j * x_j * z_j
            PVHSS_Mult(tmp, b, pk, ek, tmp, ct_x[i], prf_key);// x_i * z_i * x_j * z_j * x_j * z_j
            // appear 3 times s.a. 221 212 122
            PVHSS_cMult(tmp, b, pk, ek, tmp, 3, prf_key); // 3 * x_i * z_i * x_j * z_j * x_j * z_j
            PVHSS_ADD(term4, b, pk, ek, term4, tmp, prf_key); // term4 += 3 * x_i * z_i * x_j * z_j * x_j * z_j
        }
    } // \sum x_i * z_i * x_j * z_j * x_k * z_k when i=j or i=k or j=k
    for (auto i = 0; i < features; i++){
        for (auto j = i + 1; j < features; j++){
            for (auto k = j + 1; k < features; k++){
                if ((reg_xz_xz[i][j].s == 0 && reg_xz_xz[i][j].s_ == 0) ||  ct_x[k] == 0 || ct_z[k] == 0){
                    continue;
                }
                REG tmp;
                PVHSS_Mult(tmp, b, pk, ek, reg_xz_xz[i][j], ct_z[k], prf_key);// z_k * x_j * z_j * x_i * z_i
                PVHSS_Mult(tmp, b, pk, ek, tmp, ct_x[k], prf_key);// x_k * z_k * x_j * z_j * x_i * z_i
                // appear 6 times s.a. 123 132 213 231 312 321
                PVHSS_cMult(tmp, b, pk, ek, tmp, 6, prf_key); // 6 * x_j * z_j * x_i * z_i * x_i * z_i
                PVHSS_ADD(term4, b, pk, ek, term4, tmp, prf_key); // term4 += 6 * x_j * z_j * x_i * z_i * x_i * z_i
            }
        }
    }// \sum x_i * z_i * x_j * z_j * x_k * z_k
    PVHSS_Mult(term4, b, pk, ek, term4, ct_ggg, prf_key); // gamma^3  * (\sum x_i * z_i * x_j * z_j)

    // Add all terms together
    PVHSS_ADD(reg, b, pk, ek, reg, term1, prf_key);
    PVHSS_ADD(reg, b, pk, ek, reg, term2, prf_key);
    PVHSS_ADD(reg, b, pk, ek, reg, term3, prf_key);
    PVHSS_ADD(reg, b, pk, ek, reg, term4, prf_key);

    time = GetTime() - time;
    std::cout << "Online computing time: " << time * 1000 << " ms" << std::endl;
    bench_time << time * 1000 << std::endl;
}


// compute SVM predict
// y = \sum alpha_j * K(x_j, z) + b
void predict(REG &reg, int b, PVHSSPK pk, PVHSSEK ek, int &prf_key,
             Vec<Vec<ZZ>> ct_X, Vec<ZZ> ct_z, Vec<ZZ> ct_alpha,
             ZZ ct_ccc, ZZ ct_3ccg, ZZ ct_3cgg, ZZ ct_ggg, ZZ ct_b) {
    reg.s = 0;
    reg.s_ = 0;
    REG tmp;
    std::ofstream bench_time ("../benchmark/poly_time_" + std::to_string(ct_z.length()) + "_" + std::to_string(b) + ".txt");
    auto SV = ct_alpha.length();
    for (auto i = 0; i < SV; i++){
        std::cout << "SV: " << i << std::endl;
        if (ct_alpha[i] == 0){
            continue;
        }
        poly_kernel(tmp, b, pk, ek, ct_X[i], ct_ccc, ct_3ccg, ct_3cgg, ct_ggg, ct_z, prf_key, bench_time); // tmp = K(x, z)
        PVHSS_Mult(tmp, b, pk, ek, tmp, ct_alpha[i], prf_key); // tmp = alpha_j * K(x_j, z)
        PVHSS_ADD(reg, b, pk, ek, reg, tmp, prf_key); // reg += alpha_j * K(x_j, z)
    }
    tmp.s = 0;
    tmp.s_ = 0;
    PVHSS_Load(tmp, b, pk, ek, ct_b, prf_key);
    PVHSS_ADD(reg, b, pk, ek, reg, tmp, prf_key); // reg = \sum alpha_j * K(x_j, z) + b
    bench_time.close();
}


void get_user_inputs(Vec<Vec<ZZ>> &X, Vec<ZZ> &y, std::string in_file) {
    std::vector<std::vector<std::string>> input_data;

    read_csv(input_data, in_file);
    std::vector<std::vector<double>> f_data;
    data_preprocess(f_data, input_data);
    auto scale = pow(10, 15);
    Vec<Vec<ZZ>> ZZ_data;
    f_data2ZZ(ZZ_data, f_data, scale);

    for (auto r : ZZ_data){
        y.append(r[0]);
        Vec<ZZ> tmp;
        for (auto i = 1; i < r.length(); i++){
            tmp.append(r[i]);
        }
        X.append(tmp);
    }
}


void eval_svm(std::string in_file, std::string para_file, double gamma, double b, double c){
    int prf_key = 1;
    //Phase 1: Get all inputs
    ModelPara modelPara;
    set_model_paras(modelPara, gamma, b, c, para_file);

    std::cout << "SV: " << modelPara.SV << std::endl;
    std::cout << "features: " << modelPara.features << std::endl;
    std::cout << "gamma: " << modelPara.gamma << std::endl;
    std::cout << "b: " << modelPara.b << std::endl;
    std::cout << "c: " << modelPara.c << std::endl;

    Vec<Vec<ZZ>> Z_test;
    Vec<ZZ> y_test;
    get_user_inputs(Z_test, y_test, in_file);

    Vec<Vec<ZZ>> ct_X;
    ct_X.SetLength(modelPara.SV);
    for (auto i = 0; i < modelPara.SV; i++){
        ct_X[i].SetLength(modelPara.features);
    }
    Vec<Vec<ZZ>> ct_Z;
    ct_Z.SetLength(Z_test.length());
    for (auto i = 0; i < Z_test.length(); i++){
        ct_Z[i].SetLength(Z_test[0].length());
    }
    Vec<ZZ> ct_alpha;
    ct_alpha.SetLength(modelPara.SV);
    ZZ ct_c;
    ZZ ct_gamma;
    ZZ ct_b;



    //Phase 2: Share (HSS_Enc) all inputs
    PVHSSPK pk;
    PVHSSEK ek0, ek1;
    PVHSSPVK pvk;

    PVHSS_Gen(pk, ek0, ek1, pvk);
//    read_keys_from_file(pk, ek0, ek1, pvk);

//    read_para_from_file(ct_X, ct_alpha, ct_c, ct_gamma, ct_b);

    ZZ ct_ccc, ct_3ccg, ct_3cgg, ct_ggg;
    PVHSS_Enc(ct_ccc, pk, modelPara.c * modelPara.c * modelPara.c);
    PVHSS_Enc(ct_3ccg, pk, 3 * modelPara.c * modelPara.c * modelPara.gamma);
    PVHSS_Enc(ct_3cgg, pk, 3 * modelPara.c *modelPara.gamma * modelPara.gamma);
    PVHSS_Enc(ct_ggg, pk, modelPara.gamma * modelPara.gamma * modelPara.gamma);
//    std::cout << ct_b << std::endl;
//    std::cout << ct_c << std::endl;

    // support vectors
    std::ofstream sv_file("../data/ct_sv_"+ std::to_string(modelPara.features) +".txt");
    if (!sv_file){
        std::cerr << "Cannot open the ct_sv_"+ std::to_string(modelPara.features) +".txt file." << std::endl;
        exit(1);
    }
    for (auto i = 0; i < modelPara.SV; i++){
        std::cout << "Enc " << i << "th sv" << std::endl;
        for (auto j = 0; j < modelPara.features; j++){
//            std::cout << i << " " << j << std::endl;
            if (modelPara.sv[i][j] == 0 ){
                ct_X[i][j] = 0;
                sv_file << ct_X[i][j] << " ";
                continue;
            }
            PVHSS_Enc(ct_X[i][j], pk, modelPara.sv[i][j]);
            sv_file << ct_X[i][j] << " ";
        }
        sv_file << std::endl;
    }
    sv_file.close();

    // inputs
    std::ofstream input_time("../benchmark/input_time_dp_"+ std::to_string(modelPara.features) +".txt");
    if (!input_time) {
        std::cerr << "Cannot open the input_time_"+ std::to_string(modelPara.features) +".txt file." << std::endl;
        exit(1);
    }
    for (auto i = 0; i < Z_test.length(); i++){
        auto time = GetTime();
        for (auto j = 0; j < Z_test[0].length(); j++){
            if (Z_test[i][j] == 0 ){
                ct_Z[i][j] = 0;
                continue;
            }
            PVHSS_Enc(ct_Z[i][j], pk, Z_test[i][j]);
        }
        time = GetTime() - time;
        input_time << time * 1000 << std::endl;
    }
    input_time.close();

    // alpha_j * y_j
    std::ofstream alpha_file("../data/ct_alpha_"+ std::to_string(modelPara.features) +".txt");
    if (!alpha_file){
        std::cerr << "Cannot open the ct_alpha.txt file." << std::endl;
        exit(1);
    }
    for (auto i = 0; i < modelPara.SV; i++){
        if (modelPara.alpha[i] == 0 ){
            ct_alpha[i] = 0;
            alpha_file << ct_alpha[i] << std::endl;
            continue;
        }
        PVHSS_Enc(ct_alpha[i], pk, modelPara.alpha[i]);
        alpha_file << ct_alpha[i] << " ";
    }
    alpha_file << std::endl;
    PVHSS_Enc(ct_c, pk, modelPara.c);
    alpha_file << ct_c << std::endl;
    PVHSS_Enc(ct_gamma, pk, modelPara.gamma);
    alpha_file << ct_gamma << std::endl;
    PVHSS_Enc(ct_b, pk, modelPara.b);
    alpha_file << ct_b << std::endl;
    alpha_file.close();
//


    Vec<ZZ> V;
    V.SetLength(modelPara.SV);

/*    std::ofstream bench_time ("../benchmark/poly_time_test.txt");
    for (auto i = 0; i < modelPara.SV; i++){
        dot_prod_plaintext(V[i], modelPara.sv[i], Z_test[0]);
        power(V[i], modelPara.gamma * V[i] + modelPara.c, 3);
        std::cout << "True:" << V[i] <<std::endl;
        REG tmp0, tmp1;
        poly_kernel(tmp0, 0, pk, ek0, ct_X[i], ct_ccc, ct_3ccg, ct_3cgg, ct_ggg, ct_Z[0], prf_key, bench_time);
        poly_kernel(tmp1, 1, pk, ek1, ct_X[i], ct_ccc, ct_3ccg, ct_3cgg, ct_ggg, ct_Z[0], prf_key, bench_time);

        ZZ v, v0, v1;
        ep_t g0_tau;
        ep2_t g1_tau;
        PVHSS_Output(v0, g0_tau, g1_tau, 0, pk, ek0, tmp0, prf_key);
        PVHSS_Output(v1, g0_tau, g1_tau, 1, pk, ek1, tmp1, prf_key);

        auto pass = PVHSS_Ver(v, v0, v1, g0_tau, g1_tau, pk, pvk);
        if (!pass){
            return;
        }
    }
    bench_time.close();*/


/*
    REG tmp0, tmp1;
    poly_kernel(tmp0, 0, pk, ek0, ct_X[3], ct_ccc, ct_3ccg, ct_3cgg, ct_ggg, ct_Z[0], prf_key, bench_time);
    poly_kernel(tmp1, 1, pk, ek1, ct_X[3], ct_ccc, ct_3ccg, ct_3cgg, ct_ggg, ct_Z[0], prf_key, bench_time);

    ZZ v, v0, v1;
    ep_t g0_tau;
    ep2_t g1_tau;
    PVHSS_Output(v0, g0_tau, g1_tau, 0, pk, ek0, tmp0, prf_key);
    PVHSS_Output(v1, g0_tau, g1_tau, 1, pk, ek1, tmp1, prf_key);

    auto pass = PVHSS_Ver(v, v0, v1, g0_tau, g1_tau, pk, pvk);
    if (pass){
        std::cout << v << std::endl;
    }*/


//    // compute all dot product for some z
    Vec<ZZ> ct_V;

    ct_V.SetLength(modelPara.SV);

    std::ofstream ver_time("../benchmark/ver_time_"+ std::to_string(modelPara.features) +".txt");
    std::ofstream dot_prod_time_0("../benchmark/dot_prod_time_0_"+ std::to_string(modelPara.features) +".txt");
    std::ofstream dot_prod_time_1("../benchmark/dot_prod_time_1_"+ std::to_string(modelPara.features) +".txt");
    // std::ofstream squared_euclidean_distance_time_0("../benchmark/squared_euclidean_distance_time_0_"+ std::to_string(modelPara.features) +".txt");
    // std::ofstream squared_euclidean_distance_time_1("../benchmark/squared_euclidean_distance_time_1_"+ std::to_string(modelPara.features) +".txt");


    if (!ver_time || !squared_euclidean_distance_time_0 || !squared_euclidean_distance_time_1){
        std::cerr << "Cannot open the file." << std::endl;
        exit(1);
    }
    for (auto i = 0; i < modelPara.SV; i++){
//        std::cout << "start " << i << std::endl;
        REG tmp0, tmp1;
        squared_euclidean_distance(tmp0, 0, pk, ek0, prf_key, ct_X[i], ct_Z[0], ct_Z[0], squared_euclidean_distance_time_0);
        squared_euclidean_distance(tmp1, 1, pk, ek1, prf_key, ct_X[i], ct_Z[0], ct_Z[0], squared_euclidean_distance_time_1);
//        dot_prod(tmp0, 0, pk, ek0, prf_key, ct_X[i], ct_Z[0], dot_prod_time_0);
//        dot_prod(tmp1, 1, pk, ek1, prf_key, ct_X[i], ct_Z[0], dot_prod_time_1);
//        std::cout << "end" << std::endl;
        ZZ v, v0, v1;
        ZZ g_tau0, g_tau1;
        PVHSS_Output(v0, g_tau0, 0, pk, ek0, tmp0, prf_key);
        PVHSS_Output(v1, g_tau1, 1, pk, ek1, tmp1, prf_key);


        auto time = GetTime();
        auto pass = PVHSS_Ver(v, v0, v1, g_tau0, g_tau1, pk, pvk);
        time = GetTime() - time;
        ver_time << time * 1000 << std::endl;
        if (!pass){
            exit(1);
        }
        // Send v to user to compute poly kernel
        power(v, modelPara.gamma * v + modelPara.c, 3);
//        std::cout << "After power: " << v << std::endl;
        PVHSS_Enc(ct_V[i], pk, v);
//        dot_prod_plaintext(V[i], modelPara.sv[i], Z_test[0]);
//        power(V[i], modelPara.gamma * V[i] + modelPara.c, 3);
//        std::cout << "Plaintext after power: " << V[i] << std::endl;
    }


    std::cout << "start" << std::endl;
    REG tmp0, tmp1;
//    dot_prod(tmp0, 0, pk, ek0, prf_key, ct_V, ct_alpha, dot_prod_time_0);
//    dot_prod(tmp1, 1, pk, ek1, prf_key, ct_V, ct_alpha, dot_prod_time_1);
    // add b
    REG tmp0_b, tmp1_b;
    PVHSS_Load(tmp0_b, 0, pk, ek0, ct_b, prf_key);
    PVHSS_Load(tmp1_b, 1, pk, ek1, ct_b, prf_key);
    PVHSS_ADD(tmp0, 0, pk, ek0, tmp0, tmp0_b, prf_key);
    PVHSS_ADD(tmp1, 1, pk, ek1, tmp1, tmp1_b, prf_key);
    std::cout << "end" << std::endl;
    ZZ v, v0, v1;
    ZZ g_tau0, g_tau1;
    PVHSS_Output(v0, g_tau0, 0, pk, ek0, tmp0, prf_key);
    PVHSS_Output(v1, g_tau1, 1, pk, ek1, tmp1, prf_key);

    auto time = GetTime();
    auto pass = PVHSS_Ver(v, v0, v1, g_tau0, g_tau1, pk, pvk);
    time = GetTime() - time;
    ver_time << time * 1000 << std::endl;

    ver_time.close();
//    dot_prod_time_0.close();
//    dot_prod_time_1.close();
    if (!pass){
//        exit(1);
    }

    ZZ res;
    dot_prod_plaintext(res, V, modelPara.alpha);
    res += modelPara.b;
    std::cout << res << std::endl;
}

void dot_prod(REG &reg, int b, PVHSSPK pk, PVHSSEK ek, int &prf_key, Vec<ZZ> ct_x, Vec<ZZ> ct_z, std::ofstream &dot_prod_time) {
    Vec<REG> reg_x;
    reg_x.SetLength(ct_x.length());
    for (auto i = 0; i < ct_x.length(); i++) {
        if (ct_x[i] == 0){
            reg_x[i].s = 0;
            reg_x[i].s_ = 0;
            continue;
        }
        PVHSS_Load(reg_x[i], b, pk, ek, ct_x[i], prf_key);
    }
    REG tmp;
    std::cout << "Start online computing" << std::endl;
    double time = GetTime();
    for (auto i = 0; i < ct_x.length(); i++){
        if (ct_x[i] == 0 || ct_z[i] == 0){
            continue;
        }
        PVHSS_Mult(tmp, b, pk, ek, reg_x[i], ct_z[i], prf_key);
        PVHSS_ADD(reg, b, pk, ek, reg, tmp, prf_key);
    }
    time = GetTime() - time;
    dot_prod_time << time * 1000 << std::endl;
    std::cout << "Online computing time: " << time * 1000 << " ms" << std::endl;
}

void squared_euclidean_distance(REG &reg, int b, PVHSSPK pk, PVHSSEK ek, int &prf_key, Vec<ZZ> ct_x, Vec<ZZ> ct_z, Vec<ZZ> ct_zz, std::ofstream &squared_euclidean_distance_time){
    Vec<REG> reg_x;
    reg_x.SetLength(ct_x.length());
    for (auto i = 0; i < ct_x.length(); i++) {
        if (ct_x[i] == 0){
            reg_x[i].s = 0;
            reg_x[i].s_ = 0;
            continue;
        }
        PVHSS_Load(reg_x[i], b, pk, ek, ct_x[i], prf_key);
    }

    REG reg_xx;
    reg_xx.s = 0;
    reg_xx.s_ = 0;
    for (auto i = 0; i < ct_x.length(); i++){
        if (ct_x[i] == 0){
            continue;
        }
        REG tmp;
        PVHSS_Mult(tmp, b, pk, ek, reg_x[i], ct_x[i], prf_key);
        PVHSS_ADD(reg_xx, b, pk, ek, reg_xx, tmp, prf_key);
    }

    // Start online computing
    double time = GetTime();
    REG reg_xz;
    reg_xz.s = 0;
    reg_xz.s_ = 0;
    for (auto i = 0; i < ct_x.length(); i++){
        if (ct_x[i] == 0 || ct_z[i] == 0){
            continue;
        }
        REG tmp;
        PVHSS_Mult(tmp, b, pk, ek, reg_x[i], ct_z[i], prf_key);
        PVHSS_ADD(reg_xz, b, pk, ek, reg_xz, tmp, prf_key);
    }

    Vec<REG> reg_z;
    reg_z.SetLength(ct_z.length());
    for (auto i = 0; i < ct_z.length(); i++) {
        if (ct_z[i] == 0){
            reg_z[i].s = 0;
            reg_z[i].s_ = 0;
            continue;
        }
        PVHSS_Load(reg_z[i], b, pk, ek, ct_z[i], prf_key);
    }

    REG reg_zz;
    reg_zz.s = 0;
    reg_zz.s_ = 0;
    for (auto i = 0; i < ct_z.length(); i++){
        if (ct_zz[i] == 0){
            continue;
        }
        PVHSS_Load(reg_zz, b, pk, ek, ct_zz[i], prf_key);
//        REG tmp;
//        PVHSS_Mult(tmp, b, pk, ek, reg_z[i], ct_z[i], prf_key);
//        PVHSS_ADD(reg_zz, b, pk, ek, reg_zz, tmp, prf_key);
    }

    PVHSS_ADD(reg, b, pk, ek, reg_xx, reg_zz, prf_key);
    reg_xz.s = -2 * reg_xz.s;
    reg_xz.s_ = -2 * reg_xz.s_;
    PVHSS_ADD(reg, b, pk, ek, reg, reg_xz, prf_key);

    time = GetTime() - time;
    squared_euclidean_distance_time << time * 1000 << std::endl;
    std::cout << "Online computing time: " << time * 1000 << " ms" << std::endl;
}


void test_input_dp_time(std::string in_file, std::string type, int features){
    //Phase 1: Get all inputs
    Vec<Vec<ZZ>> Z_test;
    Vec<ZZ> y_test;
    get_user_inputs(Z_test, y_test, in_file);

    Vec<Vec<ZZ>> ct_Z;
    ct_Z.SetLength(Z_test.length());
    for (auto i = 0; i < Z_test.length(); i++){
        ct_Z[i].SetLength(Z_test[0].length());
    }


    //Phase 2: Share (HSS_Enc) all inputs
    PVHSSPK pk;
    PVHSSEK ek0, ek1;
    PVHSSPVK pvk;

    PVHSS_Gen(pk, ek0, ek1, pvk);


    // inputs
    std::ofstream input_time("../benchmark/input_time_"+ type +"_"+ std::to_string(features) +".txt");
    if (!input_time) {
        std::cerr << "Cannot open the input_time_"+ std::to_string(features) +".txt file." << std::endl;
        exit(1);
    }
    for (auto i = 0; i < Z_test.length(); i++){
        auto time = GetTime();
        for (auto j = 0; j < Z_test[0].length(); j++){
            if (Z_test[i][j] == 0 ){
                ct_Z[i][j] = 0;
                continue;
            }
            PVHSS_Enc(ct_Z[i][j], pk, Z_test[i][j]);
        }
        time = GetTime() - time;
        input_time << time * 1000 << std::endl;
    }
    input_time.close();
}

void test_input_sv_time(std::string para_file, std::string type, int features){
    //Phase 1: Get all inputs
    ModelPara modelPara;
    set_model_paras(modelPara, 1, 1, 1, para_file);



    Vec<Vec<ZZ>> ct_X;
    ct_X.SetLength(modelPara.SV);
    for (auto i = 0; i < modelPara.SV; i++){
        ct_X[i].SetLength(modelPara.features);
    }

    Vec<ZZ> ct_alpha;
    ct_alpha.SetLength(modelPara.SV);
    ZZ ct_c;
    ZZ ct_gamma;
    ZZ ct_b;



    //Phase 2: Share (HSS_Enc) all inputs
    PVHSSPK pk;
    PVHSSEK ek0, ek1;
    PVHSSPVK pvk;

    PVHSS_Gen(pk, ek0, ek1, pvk);


    // support vectors
    std::ofstream input_sv_time("../benchmark/input_time_sv_"+type+"_"+ std::to_string(modelPara.features) +".txt");
    if (!input_sv_time) {
        std::cerr << "Cannot open the input_time_"+ std::to_string(modelPara.features) +".txt file." << std::endl;
        exit(1);
    }
    for (auto i = 0; i < modelPara.SV; i++){
        std::cout << "Enc " << i << "th sv" << std::endl;
        auto time = GetTime();
        for (auto j = 0; j < modelPara.features; j++){
            if (modelPara.sv[i][j] == 0 ){
                ct_X[i][j] = 0;
                continue;
            }
            PVHSS_Enc(ct_X[i][j], pk, modelPara.sv[i][j]);
        }
        time = GetTime() - time;
        input_sv_time << time * 1000 << std::endl;
    }
}

void test_batch_verify(int size){
    //Phase 2: Share (HSS_Enc) all inputs
    PVHSSPK pk;
    PVHSSEK ek0, ek1;
    PVHSSPVK pvk;

    PVHSS_Gen(pk, ek0, ek1, pvk);

    std::ofstream verify_time("../benchmark/verify_time_size_"+ std::to_string(size) + ".txt");
    for (auto k = 0; k < 200; k++){
        ZZ V, v, v0, v1;
        ZZ Tau, tau, tau0, tau1;
        RandomBnd(v0, pk.DJpk1.N);
        RandomBnd(v1, pk.DJpk1.N);
        RandomBnd(tau1, pk.DJpk1.N);
        RandomBnd(tau0, pk.DJpk1.N);

        auto time = GetTime();
        for (auto i = 0; i < size; i++){ // number of intervals
            v = v1- v0;
            tau = tau1 - tau0;
            V += v;
            Tau += tau;
        }


        ZZ left, right;
        left = PowerMod(pk.g, Tau, pk.r_);
        right = PowerMod(pvk, V, pk.r_);
        if (left != right){
            std::cout << "Fail" << std::endl;
        } else {
            std::cout << "Pass" << std::endl;
        }
        time = GetTime() - time;
        verify_time << time * 1000 << std::endl;
    }

}

void test_verify(int size){
    PVHSSPK pk;
    PVHSSEK ek0, ek1;
    PVHSSPVK pvk;

    PVHSS_Gen(pk, ek0, ek1, pvk);

    std::ofstream verify_time("../benchmark/verify_time_poly_size_"+ std::to_string(size) + ".txt");
    for (auto i = 0; i < 200; i++){
        ZZ v, v0, v1;
        ZZ tau, tau0, tau1;
        RandomBnd(v0, pk.DJpk1.N);
        RandomBnd(v1, pk.DJpk1.N);
        RandomBnd(tau1, pk.DJpk1.N);
        RandomBnd(tau0, pk.DJpk1.N);

        auto time = GetTime();
        v = v1 - v0;
        tau = tau1 - tau0;

        ZZ left, right;
        left = PowerMod(pk.g, tau, pk.r_);
        right = PowerMod(pvk, v, pk.r_);
        if (left != right){
            std::cout << "Fail" << std::endl;
        } else {
            std::cout << "Pass" << std::endl;
        }
        time = GetTime() - time;
        verify_time << time * 1000 << std::endl;
    }
}

void test_input_size(){
    //Phase 1: Get all inputs
    Vec<Vec<ZZ>> Z_test;
    Vec<ZZ> y_test;
    get_user_inputs(Z_test, y_test, "../data/SMS_test_data_1000.csv");

    auto count = 0;
    for (auto i = 0; i < Z_test.length(); i++){
        for (auto j = 0; j < Z_test[0].length(); j++){
            if (Z_test[i][j] != 0){
                count ++;
            }
        }

        if ((i + 1)%200 == 0){
            std::cout << "Size: " << i + 1 << " -- " << count << std::endl;
        }
    }

    PVHSSPK pk;
    PVHSSEK ek0, ek1;
    PVHSSPVK pvk;

    PVHSS_Gen(pk, ek0, ek1, pvk);
}

