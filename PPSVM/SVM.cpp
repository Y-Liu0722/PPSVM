#include "SVM.h"

vec_ZZ_p HomMVMult(int b, const EK &ek, const PKE_Para &pkePara, const Vec<MemoryV> &t_M, const Ciphertext &C_v, int m, int n) {
    vec_ZZ_p u;
    MemoryV t_P;
    Vec<vec_ZZ_p> p_b;
    vec_ZZ_p p;
    u.SetLength(m);
    p_b.SetLength(t_M.length());
    for (auto i = 0; i < t_M.length(); i++){
        BKS_Mult(t_P, b, ek.bksEk, pkePara, t_M[i], C_v);
        {
            // switch to mod r
            ZZ_pPush push(pkePara.r_context);
//            std::cout << BKS_Output(b, ek.bksEk, t_P, pkePara) << std::endl;
            p_b[i] = Decode(pkePara, BKS_Output(b, ek.bksEk, t_P, pkePara));
//            std::cout << "Decode time: " << time << std::endl;
            p.append(p_b[i]);
        }
    }

    ZZ_pPush push(pkePara.r_context);
    for (auto i = 0; i < m; i++){
        for (auto j = 0; j < n; j++){
            u[i] += p[i * n + j];
        }
    }

    return u;
}

void SVM_Gen(PubPara &para, PKE_PK &pk, EK &ek1, EK &ek2, PVK &pvk) {
    BKS_Gen(ek1.bksEk, ek2.bksEk, pk, para.pkePara);

//    std::cout << "Complete the BKS Key Gen." << std::endl;

    ZZ r;
    ZZ r_;

//    GenGermainPrime(r, 3072);
    r = conv<ZZ>("42865672803402392469208087890080561026811317745152730289594047845462967569415"
                 "7674828456008147789655391216568967784547637513648556483025116956674453690415901"
                 "0911518880273586983896001426533858413249502523635626303820137243885612640280451"
                 "0057586615854068217862577095756234056971371545292003983646233796637548571023408"
                 "4645611467640072786960231757938750507620841641951901599903058088523210916434047"
                 "2652939331764346519853889743643993254625849578145672701128473372638831592592938"
                 "4556584468336489462474598221704840124379808186285103767340511905891855355679953"
                 "0640194732355189059074574407114729457247875683928981084217160937060021404284943"
                 "8394509474331176473722225505793552075745858120223098556873263659900797013766243"
                 "1684621646979901821356760154000986324685287980290257408756633297197123421381720"
                 "3589497021548742986081673903491684278872270686278621918700762494954251941588367"
                 "1947668542966362782220332615444132308131148146979535742191");
    r_ = 2*r + 1;
//
//
//    std::cout << r << std::endl;
//    std::cout << r_ << std::endl;

    ZZ g;
//    ZZ h;
//    do {
//        RandomBnd(h, r_);
//        PowerMod(g, h, 2, r_);
//    } while (g == 1);
//    std::cout << g << std::endl;
//    std::cout << h << std::endl;
//    ZZ t;
//    PowerMod(t, g, r, r_);
//    std::cout << t << std::endl;

    g = conv<ZZ>("19058280161271858772252962237914698741920826123299786961480914111132670032494"
                 "2885887018443151672176902040029652415475246969022490880264854766001262469438137"
                 "8665958625810439324805043625732198261172044942926241746519761587203855908447486"
                 "4613572176294694277726295707625133046769435392064112708332139069604125360902373"
                 "6211778742727244545621001126418227523107147563419090342230175775734820046306889"
                 "6401530503958447146581467785927621373689983779232783667462302546004317014087343"
                 "5017067672914648002383688844581367081584205898216514684965181440079127763624345"
                 "2799930518503071821874869016880898381110185534075331957205611706919659091219219"
                 "1745310947729250858122640415414402910659805357943738128654318269399551264735108"
                 "6230617333206277686781615854581769667622175363140346718592686700114671844183692"
                 "7087787862490725466548841234399362625736644717061783984965055190085746688028200"
                 "6795939730690931760052721568174553762330444231581320068924");

    para.g = g;
    para.r = r;
    para.r_ = r_;

//    std::cout << "Change mod" << std::endl;

    ZZ_p::init(para.pkePara.r);
    ZZ_p alpha;
    random(alpha);
    PowerMod(pvk, para.g, rep(alpha), para.r_);

    vec_ZZ_p vec_alpha;
    vec_alpha.SetLength(para.pkePara.N, alpha);
//    for (auto i = 0; i < para.pkePara.N; i++){
//        conv(vec_alpha[i], alpha);
//    }
    ZZ_p::init(para.pkePara.q);

    BKS_Enc(ek1.C_delta, para.pkePara, pk, Encode( para.pkePara, vec_alpha));
    ek2.C_delta = ek1.C_delta;
}

void SVM_ModelEnc_basic(Vec<Ciphertext> &C, Ciphertext &C_d, Ciphertext &C_b, PubPara &para, const PKE_PK &pk,
                        const ModelPara &modelPara) {
    C.SetLength(modelPara.n_);
    BKS_Enc(C_d, para.pkePara, pk, modelPara.hat_alpha);
    BKS_Enc(C_b, para.pkePara, pk, modelPara.hat_b);

    for (auto i = 0; i < modelPara.n_; i++) {
        BKS_Enc(C[i], para.pkePara, pk, modelPara.hat_X[i]);
    }
    para.m = modelPara.m;
    para.n = modelPara.n;
}

void SVM_ModelEnc_improved(Vec<Ciphertext> &C, Ciphertext &C_d, Ciphertext &C_b, Ciphertext &C_g, Ciphertext &C_c,
                            PubPara &para, const PKE_PK &pk, const ModelPara &modelPara) {
    C.SetLength(modelPara.n_);
    BKS_Enc(C_d, para.pkePara, pk, modelPara.hat_alpha);
    BKS_Enc(C_b, para.pkePara, pk, modelPara.hat_b);
    BKS_Enc(C_g, para.pkePara, pk, modelPara.hat_gamma);
    BKS_Enc(C_c, para.pkePara, pk, modelPara.hat_c);

    for (auto i = 0; i < modelPara.n_; i++) {
        BKS_Enc(C[i], para.pkePara, pk, modelPara.hat_X[i]);
    }
    para.m = modelPara.m;
    para.n = modelPara.n;
}

// compute protocol of basic scheme
void Compute_basic_poly(ZZ_p &y_1, ZZ_p &y_2, ZZ &g_phi_1, ZZ &g_phi_2, const EK &ek1, const EK &ek2, const PubPara &para,
                    const Vec<Ciphertext> &C, const Ciphertext &C_d, const Ciphertext &C_b, const ZZ_p &gamma,
                    const ZZ_p &c, const ZZ &p, const PKE_PK &pk, const PVK &pvk, const vec_ZZ_p &z,
                    std::vector<double> &Time) { // user inputs
    // Offline
    Vec<MemoryV> t_X_1, t_X_2;
    Vec<MemoryV> t_X_delta_1, t_X_delta_2;
    MemoryV t_b_1, t_b_2, t_d_1, t_d_2;
    t_X_1.SetLength(C.length());
    t_X_delta_1.SetLength(C.length());
    t_X_2.SetLength(C.length());
    t_X_delta_2.SetLength(C.length());
    // Server 1
    auto time = GetTime();
    for (auto i = 0; i < C.length(); i++){
        BKS_Load(t_X_1[i], 1, ek1.bksEk, para.pkePara, C[i]);
        BKS_Mult(t_X_delta_1[i], 1, ek1.bksEk, para.pkePara, t_X_1[i], ek1.C_delta);
    }
    BKS_Load(t_d_1, 1, ek1.bksEk, para.pkePara, C_d);
    BKS_Load(t_b_1, 1, ek1.bksEk, para.pkePara, C_b);
    Time[0] = GetTime() - time;
//    std::cout << "Server 1 offline time: " << Time[0] << std::endl;

    // Server 2
    time = GetTime();
    for (auto i = 0; i < C.length(); i++){
        BKS_Load(t_X_2[i], 2, ek2.bksEk, para.pkePara, C[i]);
        BKS_Mult(t_X_delta_2[i], 2, ek2.bksEk, para.pkePara, t_X_2[i], ek2.C_delta);
    }
    BKS_Load(t_d_2, 2, ek2.bksEk, para.pkePara, C_d);
    BKS_Load(t_b_2, 2, ek2.bksEk, para.pkePara, C_b);
    Time[1] = GetTime() - time;
//    std::cout << "Server 2 offline time: " << Time[1] << std::endl;

    // Online
    // User input
    vec_ZZ_p z_; // length N
    for (auto i = 0 ; i < para.pkePara.N / z.length(); i++){
        z_.append(z);
    }
    auto hat_z = Encode(para.pkePara, z_);
    time = GetTime();
    Ciphertext C_z;
    BKS_Enc(C_z, para.pkePara, pk, hat_z);
    Time[4] = GetTime() - time;
//    std::cout << "User enc time: " << Time[4] << std::endl;

    // Sending C_z to server

    // Server 1 compute
    time = GetTime();
    vec_ZZ_p u_1 = HomMVMult(1, ek1, para.pkePara, t_X_1, C_z, para.m, para.n);
    vec_ZZ_p tau_1 = HomMVMult(1, ek1, para.pkePara, t_X_delta_1, C_z, para.m, para.n);
    vec_ZZ g_tau_1;
    g_tau_1.SetLength(para.m);
    for (auto i = 0; i < para.m; i++){
        PowerMod(g_tau_1[i], para.g, rep(tau_1[i]), para.r_);
    }
    Time[2] = GetTime() - time;

    // Server 2 compute
    time = GetTime();
    vec_ZZ_p u_2 = HomMVMult(2, ek2, para.pkePara, t_X_2, C_z, para.m, para.n);
    vec_ZZ_p tau_2 = HomMVMult(2, ek2, para.pkePara, t_X_delta_2, C_z, para.m, para.n);

    vec_ZZ g_tau_2;
    g_tau_2.SetLength(para.m);
    for (auto i = 0; i < para.m; i++){
        PowerMod(g_tau_2[i], para.g, rep(tau_2[i]), para.r_);
    }
    Time[3] = GetTime() - time;

    // Sending u and tau to the user
    // User check
    vec_ZZ_p k;
    {
        ZZ_pPush push(para.pkePara.r_context);
        ZZ left, right, t_right;
        right = 1;
        ZZ_p sum_u;
        for (auto i = 0; i < para.m; i++) {
            sum_u += u_1[i] + u_2[i];
            MulMod(t_right, g_tau_1[i], g_tau_2[i], para.r_);
            MulMod(right, right, t_right, para.r_);
        }
        PowerMod(left, pvk, rep(sum_u), para.r_);

        if (left != right) {
//            std::cout << "Error when checking u" << std::endl;
        }

        // User compute kernel function
        time = GetTime();
        k.SetLength(para.m);
        for (auto i = 0; i < para.m; i++) {
            k[i] = gamma * (u_1[i] + u_2[i]) + c;
            power(k[i], k[i], p);
        }
        Time[5] = GetTime() - time;
    }
    Ciphertext C_k;
    time = GetTime();
    auto hat_k = Encode(para.pkePara, k);
    BKS_Enc(C_k, para.pkePara, pk, hat_k);
    Time[5] += GetTime() - time;

    // Sending C_k to servers
    MemoryV t_kd_1, t_kd_2, t_y_1, t_y_2, t_phi_1, t_phi_2;
    // Server 1 compute
    time = GetTime();
    BKS_Mult(t_kd_1, 1, ek1.bksEk, para.pkePara, t_d_1, C_k);
    BKS_ADD1(t_y_1, 1, ek1.bksEk, t_kd_1, t_b_1);
    BKS_Mult(t_phi_1, 1, ek1.bksEk, para.pkePara, t_y_1, ek1.C_delta);

    vec_ZZ_p y_vec_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_y_1, para.pkePara));
    vec_ZZ_p phi_vec_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_phi_1, para.pkePara));
    {
        ZZ_pPush push(para.pkePara.r_context);
        for (auto i = 0; i < para.m; i++) {
            y_1 += y_vec_1[i];
            g_phi_1 += rep(phi_vec_1[i]);
        }
        PowerMod(g_phi_1, para.g, g_phi_1, para.r_);
    }
    Time[2] += GetTime() - time;
//    std::cout << "Server 1 compute time: " << Time[2] << std::endl;

    // Server 2 compute
    time = GetTime();
    BKS_Mult(t_kd_2, 2, ek2.bksEk, para.pkePara, t_d_2, C_k);
    BKS_ADD1(t_y_2, 2, ek2.bksEk, t_kd_2, t_b_2);
    BKS_Mult(t_phi_2, 2, ek2.bksEk, para.pkePara, t_y_2, ek2.C_delta);

    vec_ZZ_p y_vec_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_y_2, para.pkePara));
    vec_ZZ_p phi_vec_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_phi_2, para.pkePara));
    {
        ZZ_pPush push(para.pkePara.r_context);
        for (auto i = 0; i < para.m; i++) {
            y_2 += y_vec_2[i];
            g_phi_2 += rep(phi_vec_2[i]);
        }
        PowerMod(g_phi_2, para.g, g_phi_2, para.r_);
    }
    Time[3] += GetTime() - time;
//    std::cout << "Server 2 compute time: " << Time[3] << std::endl;
}

void Compute_basic_rbf(ZZ_p& y_1, ZZ_p& y_2, ZZ &g_phi_1, ZZ &g_phi_2, const EK &ek1, const EK &ek2, const PubPara &para,
                        const Vec<Ciphertext> &C, const Ciphertext& C_d, const Ciphertext &C_b, const ZZ_p& gamma, // server inputs
                        const PKE_PK& pk, const PVK& pvk, const vec_ZZ_p& z,
                        std::vector<double> &Time) {
    // Offline
    Vec<MemoryV> t_X_1, t_X_2;
    Vec<MemoryV> t_XX_1, t_XX_2;
    Vec<MemoryV> t_X_delta_1, t_X_delta_2;
    Vec<MemoryV> t_XX_delta_1, t_XX_delta_2;
    MemoryV t_b_1, t_b_2, t_d_1, t_d_2;
    t_X_1.SetLength(C.length());
    t_X_delta_1.SetLength(C.length());
    t_X_2.SetLength(C.length());
    t_X_delta_2.SetLength(C.length());
    t_XX_1.SetLength(C.length());
    t_XX_delta_1.SetLength(C.length());
    t_XX_2.SetLength(C.length());
    t_XX_delta_2.SetLength(C.length());

    vec_ZZ_p h_1, h_2, eta_1, eta_2;
    vec_ZZ_p H_1, H_2, Eta_1, Eta_2;
    vec_ZZ_p h_i_1, h_i_2, eta_i_1, eta_i_2;
    H_1.SetLength(para.m);
    H_2.SetLength(para.m);
    Eta_1.SetLength(para.m);
    Eta_2.SetLength(para.m);
    auto dec_1 = 0.0;
    auto dec = 0.0;
    // Server 1
    auto time = GetTime();
    for (auto i = 0; i < C.length(); i++){
        BKS_Load(t_X_1[i], 1, ek1.bksEk, para.pkePara, C[i]);
        BKS_Mult(t_X_delta_1[i], 1, ek1.bksEk, para.pkePara, t_X_1[i], ek1.C_delta);
        BKS_Mult(t_XX_1[i], 1, ek1.bksEk, para.pkePara, t_X_1[i], C[i]);
        BKS_Mult(t_XX_delta_1[i], 1, ek1.bksEk, para.pkePara, t_X_delta_1[i], C[i]);
        dec = GetTime();
        h_i_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_XX_1[i], para.pkePara));
        eta_i_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_XX_delta_1[i], para.pkePara));
        dec_1 += GetTime() - dec;
        h_1.append(h_i_1);
        eta_1.append(eta_i_1);
    }
    {
        ZZ_pPush push(para.pkePara.r_context);
        for (auto j = 0; j < para.m; j++) {
            for (auto i = 0; i < para.n; i++) {
                H_1[j] += h_1[j * para.n + i];
                Eta_1[j] += eta_1[j * para.n + i];
            }
        }
    }
    BKS_Load(t_d_1, 1, ek1.bksEk, para.pkePara, C_d);
    BKS_Load(t_b_1, 1, ek1.bksEk, para.pkePara, C_b);
    Time[0] = GetTime() - time - dec_1;


    // Server 2
    auto dec_2 = 0.0;
    time = GetTime();
    for (auto i = 0; i < C.length(); i++){
        BKS_Load(t_X_2[i], 2, ek2.bksEk, para.pkePara, C[i]);
        BKS_Mult(t_X_delta_2[i], 2, ek2.bksEk, para.pkePara, t_X_2[i], ek2.C_delta);
        BKS_Mult(t_XX_2[i], 2, ek2.bksEk, para.pkePara, t_X_2[i], C[i]);
        BKS_Mult(t_XX_delta_2[i], 2, ek2.bksEk, para.pkePara, t_X_delta_2[i], C[i]);
        dec = GetTime();
        h_i_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_XX_2[i], para.pkePara));
        eta_i_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_XX_delta_2[i], para.pkePara));
        dec_2 += GetTime() - dec;
        h_2.append(h_i_2);
        eta_2.append(eta_i_2);
    }
    for (auto j = 0; j < para.m; j++){
        for (auto i = 0; i < para.n; i++){
            H_2[j] += h_2[j * para.n + i];
            Eta_2[j] += eta_2[j * para.n + i];
        }
    }
    BKS_Load(t_d_2, 2, ek2.bksEk, para.pkePara, C_d);
    BKS_Load(t_b_2, 2, ek2.bksEk, para.pkePara, C_b);
    Time[1] = GetTime() - time - dec_2;

    // Online
    // Online
    // User input
    vec_ZZ_p z_; // length N
    for (auto i = 0 ; i < para.pkePara.N / z.length(); i++){
        z_.append(z);
    }
    auto hat_z = Encode(para.pkePara, z_);
    time = GetTime();
    Ciphertext C_z;
    BKS_Enc(C_z, para.pkePara, pk, hat_z);
    Time[4] = GetTime() - time;

    // Sending C_z to server

    // Server 1 compute
    time = GetTime();
    vec_ZZ_p u_1 = HomMVMult(1, ek1, para.pkePara, t_X_1, C_z, para.m, para.n);
    vec_ZZ_p tau_1 = HomMVMult(1, ek1, para.pkePara, t_X_delta_1, C_z, para.m, para.n);
    vec_ZZ g_tau_1;
    g_tau_1.SetLength(para.m);
    for (auto i = 0; i < para.m; i++){
        u_1[i] = H_1[i] - 2 * u_1[i];
        tau_1[i] = Eta_1[i] - 2 * tau_1[i];
        PowerMod(g_tau_1[i], para.g, rep(tau_1[i]), para.r_);
    }
    Time[2] = GetTime() - time;

    // Server 2 compute
    time = GetTime();
    vec_ZZ_p u_2 = HomMVMult(2, ek2, para.pkePara, t_X_2, C_z, para.m, para.n);
    vec_ZZ_p tau_2 = HomMVMult(2, ek2, para.pkePara, t_X_delta_2, C_z, para.m, para.n);

    vec_ZZ g_tau_2;
    g_tau_2.SetLength(para.m);
    for (auto i = 0; i < para.m; i++){
        u_2[i] = H_2[i] - 2 * u_2[i];
        tau_2[i] = Eta_2[i] - 2 * tau_2[i];
        PowerMod(g_tau_2[i], para.g, rep(tau_2[i]), para.r_);
    }
    Time[3] = GetTime() - time;

    // Sending u and tau to the user
    // User check
    vec_ZZ_p k;
    {

        ZZ_pPush push(para.pkePara.r_context);
        ZZ left, right, t_right;
        right = 1;
        ZZ_p sum_u;
        for (auto i = 0; i < para.m; i++) {
            sum_u += u_1[i] + u_2[i];
            MulMod(t_right, g_tau_1[i], g_tau_2[i], para.r_);
            MulMod(right, right, t_right, para.r_);
        }
        PowerMod(left, pvk, rep(sum_u), para.r_);

        if (left != right) {
//            std::cout << "Error when checking u" << std::endl;
        }

        // User compute kernel function
        time = GetTime();
        k.SetLength(para.m);
        ZZ_p sum_z;
        for (auto i = 0; i < para.n; i++) {
            sum_z += z[i] * z[i];
        }
        for (auto i = 0; i < para.m; i++) {
            k[i] = gamma * (u_1[i] + u_2[i] + sum_z);

            //TODO: exp need rescale
        }
        Time[5] += GetTime() - time;
    }
    Ciphertext C_k;
    auto hat_k = Encode(para.pkePara, k);
    time = GetTime();
    BKS_Enc(C_k, para.pkePara, pk, hat_k);
    Time[5] += GetTime() - time;
//    std::cout << "User compute time: " << Time[5] << std::endl;

    // Sending C_k to servers
    MemoryV t_kd_1, t_kd_2, t_y_1, t_y_2, t_phi_1, t_phi_2;
    // Server 1 compute
    time = GetTime();
    BKS_Mult(t_kd_1, 1, ek1.bksEk, para.pkePara, t_d_1, C_k);
    BKS_ADD1(t_y_1, 1, ek1.bksEk, t_kd_1, t_b_1);
    BKS_Mult(t_phi_1, 1, ek1.bksEk, para.pkePara, t_y_1, ek1.C_delta);

    vec_ZZ_p y_vec_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_y_1, para.pkePara));
    vec_ZZ_p phi_vec_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_phi_1, para.pkePara));
    {
        ZZ_pPush push(para.pkePara.r_context);
        for (auto i = 0; i < para.m; i++) {
            y_1 += y_vec_1[i];
            g_phi_1 += rep(phi_vec_1[i]);
        }
        PowerMod(g_phi_1, para.g, g_phi_1, para.r_);
    }
    Time[2] += GetTime() - time;
//    std::cout << "Server 1 compute time: " << Time[2] << std::endl;

    // Server 2 compute
    time = GetTime();
    BKS_Mult(t_kd_2, 2, ek2.bksEk, para.pkePara, t_d_2, C_k);
    BKS_ADD1(t_y_2, 2, ek2.bksEk, t_kd_2, t_b_2);
    BKS_Mult(t_phi_2, 2, ek2.bksEk, para.pkePara, t_y_2, ek2.C_delta);

    vec_ZZ_p y_vec_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_y_2, para.pkePara));
    vec_ZZ_p phi_vec_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_phi_2, para.pkePara));
    {
        ZZ_pPush push(para.pkePara.r_context);
        for (auto i = 0; i < para.m; i++) {
            y_2 += y_vec_2[i];
            g_phi_2 += rep(phi_vec_2[i]);
        }
        PowerMod(g_phi_2, para.g, g_phi_2, para.r_);
    }
    Time[3] += GetTime() - time;
}

void Compute_improved(ZZ_p& y_1, ZZ_p& y_2, ZZ &g_phi_1, ZZ &g_phi_2, const EK &ek1, const EK &ek2, const PubPara &para,
                        const Vec<Ciphertext> &C, const Ciphertext& C_d, const Ciphertext &C_b,
                        const Ciphertext& C_gamma, const Ciphertext& C_c, const ZZ& p,// server inputs
                        const PKE_PK& pk, const vec_ZZ_p& z,
                        std::vector<double> &Time){
    // Offline
    Vec<MemoryV> t_gamma_X_1, t_gamma_X_2;
    Vec<MemoryV> t_gamma_X_delta_1, t_gamma_X_delta_2;
    MemoryV t_b_1, t_b_2, t_gamma_1, t_gamma_2, t_c_1, t_c_2, t_delta_b_1, t_delta_b_2, t_delta_c_1, t_delta_c_2;
    t_gamma_X_1.SetLength(C.length());
    t_gamma_X_delta_1.SetLength(C.length());
    t_gamma_X_2.SetLength(C.length());
    t_gamma_X_delta_2.SetLength(C.length());

    // Server 1
    auto time = GetTime();
    BKS_Load(t_b_1, 1, ek1.bksEk, para.pkePara, C_b);
    BKS_Load(t_gamma_1, 1, ek1.bksEk, para.pkePara, C_gamma);
    BKS_Load(t_c_1, 1, ek1.bksEk, para.pkePara, C_c);
    BKS_Mult(t_delta_b_1, 1, ek1.bksEk, para.pkePara, t_b_1, ek1.C_delta);
    BKS_Mult(t_delta_c_1, 1, ek1.bksEk, para.pkePara, t_c_1, ek1.C_delta);
    for (auto i = 0; i < C.length(); i++){
        BKS_Mult(t_gamma_X_1[i], 1, ek1.bksEk, para.pkePara, t_gamma_1, C[i]);
        BKS_Mult(t_gamma_X_delta_1[i], 1, ek1.bksEk, para.pkePara, t_gamma_X_1[i], ek1.C_delta);
    }
    Time[0] = GetTime() - time;
    vec_ZZ_p c_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_c_1, para.pkePara));
    vec_ZZ_p theta_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_delta_c_1, para.pkePara));

    // Server 2
    time = GetTime();
    BKS_Load(t_b_2, 2, ek2.bksEk, para.pkePara, C_b);
    BKS_Load(t_gamma_2, 2, ek2.bksEk, para.pkePara, C_gamma);
    BKS_Load(t_c_2, 2, ek2.bksEk, para.pkePara, C_c);
    BKS_Mult(t_delta_b_2, 2, ek2.bksEk, para.pkePara, t_b_2, ek2.C_delta);
    BKS_Mult(t_delta_c_2, 2, ek2.bksEk, para.pkePara, t_c_2, ek2.C_delta);
    for (auto i = 0; i < C.length(); i++){
        BKS_Mult(t_gamma_X_2[i], 2, ek2.bksEk, para.pkePara, t_gamma_2, C[i]);
        BKS_Mult(t_gamma_X_delta_2[i], 2, ek2.bksEk, para.pkePara, t_gamma_X_2[i], ek2.C_delta);
    }
    Time[1] = GetTime() - time;
    vec_ZZ_p c_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_c_2, para.pkePara));
    vec_ZZ_p theta_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_delta_c_2, para.pkePara));

    // Online
    // User input
    vec_ZZ_p z_; // length N
    for (auto i = 0 ; i < para.pkePara.N / z.length(); i++){
        z_.append(z);
    }
    auto hat_z = Encode(para.pkePara, z_);
    time = GetTime();
    Ciphertext C_z;
    BKS_Enc(C_z, para.pkePara, pk, hat_z);
    Time[4] = GetTime() - time;

    // Sending C_z to server

    // Server 1 compute
    time = GetTime();
    vec_ZZ_p u_1 = HomMVMult(1, ek1, para.pkePara, t_gamma_X_1, C_z, para.m, para.n);
    vec_ZZ_p tau_1 = HomMVMult(1, ek1, para.pkePara, t_gamma_X_delta_1, C_z, para.m, para.n);
//    vec_ZZ g_tau_1;
//    g_tau_1.SetLength(para.m);
    {
        ZZ_pPush push(para.pkePara.r_context);
        for (auto i = 0; i < para.m; i++) {
            u_1[i] = u_1[i] + c_1[0];
            tau_1[i] = tau_1[i] + theta_1[0];
        }
    }

    Ciphertext C_u_1, C_tau_1;
    auto hat_u_1 = Encode(para.pkePara, u_1);
    auto hat_tau_1 = Encode(para.pkePara, tau_1);

    BKS_Enc(C_u_1, para.pkePara, pk, hat_u_1);
    BKS_Enc(C_tau_1, para.pkePara, pk, hat_tau_1);
    Time[2] = GetTime() - time;

    // Server 2 compute
    time = GetTime();
    vec_ZZ_p u_2 = HomMVMult(2, ek2, para.pkePara, t_gamma_X_2, C_z, para.m, para.n);
    vec_ZZ_p tau_2 = HomMVMult(2, ek2, para.pkePara, t_gamma_X_delta_2, C_z, para.m, para.n);

    {
        ZZ_pPush push(para.pkePara.r_context);
        for (auto i = 0; i < para.m; i++) {
            u_2[i] = u_2[i] + c_2[0];
            tau_2[i] = tau_2[i] + theta_2[0];
        }
    }

    Ciphertext C_u_2, C_tau_2;
    auto hat_u_2 = Encode(para.pkePara, u_2);
    auto hat_tau_2 = Encode(para.pkePara, tau_2);

    BKS_Enc(C_u_2, para.pkePara, pk, hat_u_2);
    BKS_Enc(C_tau_2, para.pkePara, pk, hat_tau_2);
    Time[3] = GetTime() - time;

    // Sending
    Ciphertext C_u, C_tau;
    BKS_ADD2(C_u, 1, ek1.bksEk, C_u_1, C_u_2);
    BKS_ADD2(C_tau, 1, ek1.bksEk, C_tau_1, C_tau_2);

    // Server 1
    time = GetTime();
    MemoryV t_u_1, t_u_tau_1;
    BKS_Load(t_u_1, 1, ek1.bksEk, para.pkePara, C_u);
    for (auto i = 1; i < p; i ++){
        if (i == p-1){
            BKS_Mult(t_u_tau_1, 1, ek1.bksEk, para.pkePara, t_u_1, C_tau);
        }
        BKS_Mult(t_u_1, 1, ek1.bksEk, para.pkePara, t_u_1, C_u);
    }

    MemoryV t_du_1, t_du_tau_1;
    BKS_Mult(t_du_1, 1, ek1.bksEk, para.pkePara, t_u_1, C_d);
    BKS_Mult(t_du_tau_1, 1, ek1.bksEk, para.pkePara, t_u_tau_1, C_d);

    MemoryV t_y_1, t_phi_1;
    BKS_ADD1(t_y_1, 1, ek1.bksEk, t_du_1, t_b_1);
    BKS_Mult(t_phi_1, 1, ek1.bksEk, para.pkePara, t_y_1, ek1.C_delta);

    vec_ZZ_p y_vec_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_y_1, para.pkePara));
    vec_ZZ_p phi_vec_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_phi_1, para.pkePara));
    {
        ZZ_pPush push(para.pkePara.r_context);
        for (auto i = 0; i < para.m; i++) {
            y_1 += y_vec_1[i];
            g_phi_1 += rep(phi_vec_1[i]);
        }
        PowerMod(g_phi_1, para.g, g_phi_1, para.r_);
    }
    Time[2] += GetTime() - time;

    // Server 2
    time = GetTime();
    MemoryV t_u_2, t_u_tau_2;
    BKS_Load(t_u_2, 2, ek2.bksEk, para.pkePara, C_u);
    for (auto i = 1; i < p; i ++){
        if (i == p-1){
            BKS_Mult(t_u_tau_2, 2, ek2.bksEk, para.pkePara, t_u_2, C_tau);
        }
        BKS_Mult(t_u_2, 2, ek2.bksEk, para.pkePara, t_u_2, C_u);
    }

    MemoryV t_du_2, t_du_tau_2;
    BKS_Mult(t_du_2, 2, ek2.bksEk, para.pkePara, t_u_2, C_d);
    BKS_Mult(t_du_tau_2, 2, ek2.bksEk, para.pkePara, t_u_tau_2, C_d);

    MemoryV t_y_2, t_phi_2;
    BKS_ADD1(t_y_2, 2, ek2.bksEk, t_du_2, t_b_2);
    BKS_Mult(t_phi_2, 2, ek2.bksEk, para.pkePara, t_y_2, ek2.C_delta);

    vec_ZZ_p y_vec_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_y_2, para.pkePara));
    vec_ZZ_p phi_vec_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_phi_2, para.pkePara));

    {
        ZZ_pPush push(para.pkePara.r_context);
        for (auto i = 0; i < para.m; i++) {
            y_2 += y_vec_2[i];
            g_phi_2 += rep(phi_vec_2[i]);
        }
        PowerMod(g_phi_2, para.g, g_phi_2, para.r_);
    }
    Time[3] += GetTime() - time;
}
