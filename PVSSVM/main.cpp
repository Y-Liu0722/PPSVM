#include "SDPC.h"
#include "SSEDC.h"

void test_dc_input(int features);
void test_sed_input(int features);
void test_batch_verify(int size);

int main(){



    return 0;
}

void test_sed_input(int features){
    BGN_PK pk;
    BGN_SK sk;
    bn_t **X = nullptr;
    bn_t *y = nullptr;
    size_t size;

    SSEDC_KeyGen(pk, sk);

    ModelPara modelPara;
    bn_new(modelPara.scale);
    bn_new(modelPara.gamma);
    bn_new(modelPara.b);
    bn_new(modelPara.c);
    set_model_paras(modelPara, 0.005, -1.4690208, 2, "../data/SMS_para_rbf_"+ std::to_string(features) +".csv");

    get_user_inputs(X, y, "../data/SMS_test_data_rbf_"+ std::to_string(features)+".csv", size);
    std::ofstream fgen_time("../benchmark/fgen_sed_time_"+ std::to_string(features)+".txt");
    std::ofstream sgen_time("../benchmark/sgen_sed_time_"+ std::to_string(features)+".txt");
    std::ofstream compute_time("../benchmark/compute_sed_time_"+ std::to_string(features)+".txt");

//    for (auto k = 0; k < size; k++){
//        int n = 3; // 3 server
//        int t = 1; // 1 collude most
//
//        auto shares = new bn_t*[modelPara.features + 1];
//        for (int i = 0; i < modelPara.features + 1; i++){
//            shares[i] = new bn_t[n];
//            for (int j = 0; j < n; j++){
//                bn_null(shares[i][j]);
//                bn_new(shares[i][j]);
//            }
//        }
//        bn_t *index = new bn_t[n];
//        for (int i = 0; i < n; i++){
//            bn_null(index[i]);
//            bn_new(index[i]);
//        }
//        g1_t pvk;
//        g1_new(pvk);
//
//        auto time = GetTime();
//        SSEDC_SGen(shares, index, pvk, X[k], pk, n, t, modelPara.features);
//        time = GetTime() - time;
//        sgen_time << time * 1000 << std::endl;
//    }

    for (int k = 0; k < modelPara.SV; k++){
        auto *ct_W1 = new BGN_CT[modelPara.features];
        auto *ct_W2 = new BGN_CT[modelPara.features];
        auto time = GetTime();
        SSEDC_FGen(ct_W1, ct_W2, modelPara.sv[k], pk, modelPara.features);
        time = GetTime() - time;
        fgen_time << time * 1000 << std::endl;

        int n = 3; // 3 server
        int t = 1; // 1 collude most

        auto shares = new bn_t*[modelPara.features + 1];
        for (int i = 0; i < modelPara.features + 1; i++){
            shares[i] = new bn_t[n];
            for (int j = 0; j < n; j++){
                bn_null(shares[i][j]);
                bn_new(shares[i][j]);
            }
        }
        bn_t *index = new bn_t[n];
        for (int i = 0; i < n; i++){
            bn_null(index[i]);
            bn_new(index[i]);
        }
        g1_t pvk;
        g1_new(pvk);

        time = GetTime();
        SSEDC_SGen(shares, index, pvk, X[0], pk, n, t, modelPara.features);
        time = GetTime() - time;
        sgen_time << time * 1000 << std::endl;

        // shares_rev[i] sends to server i
        auto shares_rev = new bn_t*[n];
        for (int i = 0; i < n; i++){
            shares_rev[i] = new bn_t[modelPara.features + 1];
            for (int j = 0; j < modelPara.features + 1; j++){
                bn_null(shares_rev[i][j]);
                bn_new(shares_rev[i][j]);
                bn_copy(shares_rev[i][j], shares[j][i]);
            }
        }

        auto* V = new BGN_CT[n];
        auto* O = new BGN_CT[n];
        std::cout << "***************** SV: " << k << " *****************" << std::endl;
        for (int i = 0; i < n; i ++){
            auto time = GetTime();
            SSEDC_Compute(V[i], O[i], pk, shares_rev[i], ct_W1, ct_W2, i, modelPara.features);
            time = GetTime() - time;
            compute_time << time * 1000 << std::endl;
            std::cout <<"Server " << i << ": Online computing time " << time * 1000 << " ms" << std::endl;
        }
    }
}

void test_dc_input(int features){
    BGN_PK pk;
    BGN_SK sk;
    bn_t **X = nullptr;
    bn_t *y = nullptr;
    size_t size;

    SDPC_KeyGen(pk, sk);

    ModelPara modelPara;
    bn_new(modelPara.scale);
    bn_new(modelPara.gamma);
    bn_new(modelPara.b);
    bn_new(modelPara.c);
    set_model_paras(modelPara, 0.005, -1.4690208, 2, "../data/SMS_para_"+ std::to_string(features)+".csv");



    get_user_inputs(X, y, "../data/SMS_test_data_"+ std::to_string(features)+".csv", size);
    std::ofstream fgen_time("../benchmark/fgen_dc_time_"+ std::to_string(features)+".txt");
    std::ofstream sgen_time("../benchmark/sgen_dc_time_"+ std::to_string(features)+".txt");
    std::ofstream compute_time("../benchmark/compute_dc_time_"+ std::to_string(features)+".txt");
//
//    for (auto k = 0; k < size; k++){
//        int n = 3; // 3 server
//        int t = 1; // 1 collude most
//
//        auto shares = new bn_t*[modelPara.features + 1];
//        for (int i = 0; i < modelPara.features + 1; i++){
//            shares[i] = new bn_t[n];
//            for (int j = 0; j < n; j++){
//                bn_null(shares[i][j]);
//                bn_new(shares[i][j]);
//            }
//        }
//        bn_t *index = new bn_t[n];
//        for (int i = 0; i < n; i++){
//            bn_null(index[i]);
//            bn_new(index[i]);
//        }
//        g1_t pvk;
//        g1_new(pvk);
//
//        auto time = GetTime();
//        SSEDC_SGen(shares, index, pvk, X[k], pk, n, t, modelPara.features);
//        time = GetTime() - time;
//        sgen_time << time * 1000 << std::endl;
//    }


    for (int k = 0; k < modelPara.SV; k++){
        auto *ct_W = new BGN_CT[modelPara.features];
        auto time = GetTime();
        SDPC_FGen(ct_W, modelPara.sv[k], pk, modelPara.features);
        time = GetTime() - time;
        fgen_time << time * 1000 << std::endl;

        int n = 3; // 3 server
        int t = 1; // 1 collude most

        auto shares = new bn_t*[modelPara.features + 1];
        for (int i = 0; i < modelPara.features + 1; i++){
            shares[i] = new bn_t[n];
            for (int j = 0; j < n; j++){
                bn_null(shares[i][j]);
                bn_new(shares[i][j]);
            }
        }
        bn_t *index = new bn_t[n];
        for (int i = 0; i < n; i++){
            bn_null(index[i]);
            bn_new(index[i]);
        }
        g1_t pvk;
        g1_new(pvk);

        time = GetTime();
        SDPC_SGen(shares, index, pvk, X[0], pk, n, t, modelPara.features);
        time = GetTime() - time;
        sgen_time << time * 1000 << std::endl;

        // shares_rev[i] sends to server i
        auto shares_rev = new bn_t*[n];
        for (int i = 0; i < n; i++){
            shares_rev[i] = new bn_t[modelPara.features + 1];
            for (int j = 0; j < modelPara.features + 1; j++){
                bn_null(shares_rev[i][j]);
                bn_new(shares_rev[i][j]);
                bn_copy(shares_rev[i][j], shares[j][i]);
            }
        }

        auto* V = new BGN_CT[n];
        auto* O = new BGN_CT[n];
        std::cout << "***************** SV: " << k << " *****************" << std::endl;
        for (int i = 0; i < n; i ++){
            auto time = GetTime();
            SDPC_Compute(V[i], O[i], shares_rev[i], ct_W, i, modelPara.features);
            time = GetTime() - time;
            compute_time << time * 1000 << std::endl;
            std::cout <<"Server " << i << ": Online computing time " << time * 1000 << " ms" << std::endl;
        }
    }
}

void test_batch_verify(int size){
    BGN_PK pk;
    BGN_SK sk;
    bn_t **X = nullptr;
    bn_t *y = nullptr;

    SDPC_KeyGen(pk, sk);


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
    ZZ alpha, pvk;
    RandomBnd(alpha, r);
    PowerMod(pvk, g, alpha, r_);
//    r_ = 2*r + 1;
//    bn_mul_dig(r_, r, 2);
//    bn_add_dig(r_, r_, 1);
//
//
//    do {
//        bn_rand_mod(h, r_);
//        bn_mxp_dig(g, h, 2, r_);
//    } while (bn_cmp_dig(g, 1) == RLC_EQ);
//
//    bn_t alpha;
//    bn_new(alpha);
//    bn_rand_mod(alpha, r);
//
//    bn_t pvk;
//    bn_new(pvk);
//    bn_mxp(pvk, g, alpha, r_);

    bn_t mod;
    bn_new(mod);
    bn_set_dig(mod, 524287);
    std::ofstream verify_time("../benchmark/verify_time_size_"+ std::to_string(size) +".txt");
    bn_t *index = new bn_t[3];
    for (int i = 0; i < 3; i++){
        bn_null(index[i]);
        bn_new(index[i]);
        bn_set_dig(index[i], i+1);
    }

    auto *V = new bn_t[3];
    auto *O = new bn_t[3];
    for (int i = 0; i < 3; i++){
        bn_null(V[i]);
        bn_null(O[i]);
        bn_new(V[i]);
        bn_new(O[i]);
    }

    for (int k = 0; k < 200; k++){
        auto time = GetTime();
        bn_t phi_sum, psi_sum;
        bn_new(phi_sum);
        bn_new(psi_sum);
        bn_zero(phi_sum);
        bn_zero(psi_sum);
        // verify
        for (auto i = 0; i < size; i++){
            bn_t tmp;
            bn_new(tmp);
            bn_rand_mod(tmp, mod);
            make_random_shares(V, index, tmp, mod, 3, 1);
            bn_mul_dig(tmp, tmp, 12334);
            make_random_shares(O, index, tmp, mod, 3, 1);
            for (int j = 0; j < 3; j++){
//                bn_print(V[j]);
//                bn_print(O[j]);
                Sim_Dec(V[j], V[j]);
                Sim_Dec(O[j], O[j]);
//                bn_print(V[j]);
//                bn_print(O[j]);
            }

            bn_t zero;
            bn_new(zero);
            bn_zero(zero);
            bn_t phi, psi;
            bn_new(phi);
            bn_new(psi);
            lagrange_interpolate(phi, zero, index, V, mod, 1);
            lagrange_interpolate(psi, zero, index, O, mod, 1);
//
//            std::cout << "phi: ";
//            bn_print(phi);
//            std::cout << "psi: ";
//            bn_print(psi);

            bn_rand_mod(phi, pk.q);
            bn_rand_mod(psi, pk.q);

            bn_add(phi_sum, phi_sum, phi);
            bn_add(psi_sum, psi_sum, psi);
        }

        ZZ left_side, right_side;
        ZZ phi_ZZ, psi_ZZ;
        bn2ZZ(phi_ZZ, phi_sum);
        bn2ZZ(psi_ZZ, psi_sum);
//        std::cout << "phi: " << phi_ZZ << std::endl;
//        std::cout << "psi: " << psi_ZZ << std::endl;
        PowerMod(left_side, g, phi_ZZ, r_);
        PowerMod(right_side, pvk, psi_ZZ, r_);

        if (left_side == right_side){
            std::cout << "Verify successfully" << std::endl;
        } else {
            std::cout << "Verify failed" << std::endl;
        }

        time = GetTime() - time;
        verify_time << time * 1000 << std::endl;
    }
}

void test_input_size(){
    BGN_PK pk;
    BGN_SK sk;
    bn_t **X = nullptr;
    bn_t *y = nullptr;
    size_t size;

    SDPC_KeyGen(pk, sk);

    ModelPara modelPara;
    bn_new(modelPara.scale);
    bn_new(modelPara.gamma);
    bn_new(modelPara.b);
    bn_new(modelPara.c);
    set_model_paras(modelPara, 0.005, -1.4690208, 2, "../data/SMS_para_"+ std::to_string(200)+".csv");

    get_user_inputs(X, y, "../data/SMS_test_data_"+ std::to_string(200)+".csv", size);
}

