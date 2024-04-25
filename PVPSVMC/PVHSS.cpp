#include "PVHSS.h"


void PVHSS_Gen(PVHSSPK &pk, PVHSSEK &ek0, PVHSSEK &ek1, PVHSSPVK &pvk) {
//    // init
//    bn_new(pk.g1_ord);
//    bn_new(pk.g2_ord);
//    ep_new(pk.g1_gen);
//    ep2_new(pk.g2_gen);
//    fp12_new(pk.gT_gen);
//    fp12_new(pvk);
//    fp12_new(pk.const_gt);
//    core_init();
//    ep_curve_init();
//    ep_param_set(B12_P381);
//    ep2_curve_init();
//    ep2_curve_set_twist(EP_MTYPE);
//    ep_curve_get_gen(pk.g1_gen);
//    ep2_curve_get_gen(pk.g2_gen);
//    pp_map_oatep_k12(pk.gT_gen, pk.g1_gen, pk.g2_gen);

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
    pk.g = g;
    pk.r = r;
    pk.r_ = r_;

    // compute ek
    pk.K = 1024;
    pk.s = 3;
    ZZ phi, phi_;
    DJ_Gen(pk.DJpk1, phi, pk.s);
    DJ_Gen(pk.DJpk2, phi_, pk.s-2);

    ZZ mu, nu;
    InvMod(mu, phi, pk.DJpk2.Ns);

    // N'^s' mod phi'
    ZZ N_phi;
    rem(N_phi, pk.DJpk2.Ns, phi_);
    InvMod(nu, N_phi, phi_);

    // c'
    ZZ ct_;
    DJ_Enc(ct_, pk.DJpk2, mu);
    pk.ct_ = ct_;

    // NN'2^k
    ZZ NN_, v2k, NN_2k;
    mul(NN_, pk.DJpk1.N, pk.DJpk2.N);
    power(v2k, 2, 128);
    mul(NN_2k, NN_, v2k);

    // phi_0 phi_0'
    ZZ phi0, phi0_, phi1, phi1_;
    RandomBnd(phi0, NN_2k);
    RandomBnd(phi0_, NN_2k);
//    AddMod(phi1, phi0, phi, pk.DJpk1.Ns);
//    AddMod(phi1_, phi0_, phi*nu, pk.DJpk1.Ns);
    phi1 = phi0 + phi;
    phi1_ = phi0_ + phi * nu;

    // set EKs
    ek0.phi = phi0;
    ek0.phi_ = phi0_;
//    ep_copy(ek0.g1_gen, pk.g1_gen);
    ek1.phi = phi1;
    ek1.phi_ = phi1_;
//    ep2_copy(ek1.g2_gen, pk.g2_gen);

    // compute PVK
    ZZ C_A;


    // alpha
    ZZ alpha;
    RandomBits(alpha, 1024);
//    std::cout << "alpha: " << alpha << std::endl;
    PVHSS_Enc(C_A, pk, alpha);
    ek0.C_A = C_A;
    ek1.C_A = C_A;

    PowerMod(pvk, pk.g, alpha, pk.r_);

//    power2(pk.DJpk1.Nss, 1024);
//    bn_t A_bn;
//    ep_t temp;
//    bn_new(A_bn);
//    ep_new(temp);
//    ZZ2bn(A_bn, A);
//    ep_mul_gen(temp, A_bn);
//    pp_map_oatep_k12(pvk, temp, pk.g2_gen);
//
//    ep_curve_get_ord(pk.g1_ord);
//    ep2_curve_get_ord(pk.g2_ord);
//    int size = bn_size_str(pk.g1_ord, 10);
//    char *g1_order_str = new char[size];
//    bn_write_str(g1_order_str, size, pk.g1_ord, 10);
//    ZZ g1_ord_ZZ = conv<ZZ>(g1_order_str);
//    pk.g1_ord_ZZ = g1_ord_ZZ;
//    size = bn_size_str(pk.g2_ord, 10);
//    char *g2_ord_str = new char[size];
//    bn_write_str(g2_ord_str, size, pk.g2_ord, 10);
//    ZZ g2_ord_ZZ = conv<ZZ>(g2_ord_str);
//    pk.g2_ord_ZZ = g2_ord_ZZ;
}

void PVHSS_Enc(ZZ &I, PVHSSPK pk, ZZ x) {
    DJ_Enc(I, pk.DJpk1, x);
}

void PVHSS_Load(REG &reg, int b, PVHSSPK pk, PVHSSEK ek, ZZ ct_x, int &prf_key) {
    // (C^x)^phi_b and (C^x)^phi'_b
    ZZ ct_phi, ct_phi_;
    PowerMod(ct_phi, ct_x, ek.phi, pk.DJpk1.Nss);
    PowerMod(ct_phi_, ct_x, ek.phi_, pk.DJpk1.Nss);

    // s_x^b and s'_x^b
    DJ_Dist(reg.s, pk.DJpk1, ct_phi);
    DJ_Dist(reg.s_, pk.DJpk1, ct_phi_);

    AddMod(reg.s, reg.s, PRF_ZZ(prf_key, pk.DJpk1.Ns), pk.DJpk1.Ns);
    AddMod(reg.s_, reg.s_, PRF_ZZ(prf_key, pk.DJpk1.Ns), pk.DJpk1.Ns);
}

void PVHSS_ADD(REG &reg, int b, PVHSSPK pk, PVHSSEK ek, REG reg_x, REG reg_y, int &prf_key) {
    AddMod(reg.s, reg_x.s, reg_y.s, pk.DJpk1.Ns);
    AddMod(reg.s_, reg_x.s_, reg_y.s_, pk.DJpk1.Ns);

    AddMod(reg.s, reg.s, PRF_ZZ(prf_key, pk.DJpk1.Ns), pk.DJpk1.Ns);
    AddMod(reg.s_, reg.s_, PRF_ZZ(prf_key, pk.DJpk1.Ns), pk.DJpk1.Ns);
}

void PVHSS_cMult(REG &reg, int b, PVHSSPK pk, PVHSSEK ek, REG reg_x, long long c, int &prf_key) {
    MulMod(reg.s, reg_x.s, c, pk.DJpk1.Ns);
    MulMod(reg.s_, reg_x.s_, c, pk.DJpk1.Ns);

    AddMod(reg.s, reg.s, PRF_ZZ(prf_key, pk.DJpk1.Ns), pk.DJpk1.Ns);
    AddMod(reg.s_, reg.s_, PRF_ZZ(prf_key, pk.DJpk1.Ns), pk.DJpk1.Ns);
}

void PVHSS_Mult(REG &reg, int b, PVHSSPK pk, PVHSSEK ek, REG reg_x, ZZ ct_y, int &prf_key) {
    // (C^y)^(s_b^x) and (C^y)^(s_b^x')
    ZZ ct_s, ct_s_;
    PowerMod(ct_s, ct_y, reg_x.s, pk.DJpk1.Nss);
    PowerMod(ct_s_, ct_y, reg_x.s_, pk.DJpk1.Nss);

    // s_xy^b and s'_xy^b
    DJ_Dist(reg.s, pk.DJpk1, ct_s);
    DJ_Dist(reg.s_, pk.DJpk1, ct_s_);

    AddMod(reg.s, reg.s, PRF_ZZ(prf_key, pk.DJpk1.Ns), pk.DJpk1.Ns);
    AddMod(reg.s_, reg.s_, PRF_ZZ(prf_key, pk.DJpk1.Ns), pk.DJpk1.Ns);
}

void PVHSS_Output(ZZ &v, ZZ &g_tau, int b, PVHSSPK pk, PVHSSEK ek, REG reg_x, int &prf_key) {
//    ZZ ct_s, ct_s_;
//    PowerMod(ct_s, ek.C_A, s_x[0], pk.DJpk1.Nss);
//    PowerMod(ct_s_, ek.C_A, s_x[1], pk.DJpk1.Nss);
//
//    Vec<ZZ> s;
//    s.SetLength(2);
//    // s_Ax^b and s'_Ax^b
//    DJ_Dist(s[0], pk.DJpk1, ct_s);
//    DJ_Dist(s[1], pk.DJpk1, ct_s_);
//
//    AddMod(s[0], s[0], PRF_ZZ(prf_key, pk.DJpk1.Ns), pk.DJpk1.Ns);
//    AddMod(s[1], s[1], PRF_ZZ(prf_key, pk.DJpk1.Ns), pk.DJpk1.Ns);

    REG reg;
    PVHSS_Mult(reg, b, pk, ek, reg_x, ek.C_A, prf_key);

    // v_b^x, N'^(s') s'_b^x, s_b^x - N'^(s') s'_b^x, c'^(s_b^x - N'^(s') s'_b^x )
    ZZ c_s_Ns_s_tau, tau;
    PowerMod(c_s_Ns_s_tau, pk.ct_, reg.s - pk.DJpk2.Ns * reg.s_, pk.DJpk2.Nss);
    DJ_Dist(tau, pk.DJpk2, c_s_Ns_s_tau);

    AddMod(tau, tau, PRF_ZZ(prf_key, pk.DJpk2.Ns), pk.DJpk2.Ns);

    ZZ c_s_Ns_s;
//    MulMod(Ns_s, pk.DJpk2.Ns, s_x[1], pk.DJpk1.Ns);
//    SubMod(s_Ns_s, s_x[0], Ns_s, pk.DJpk1.Ns);
    PowerMod(c_s_Ns_s, pk.ct_, reg_x.s - pk.DJpk2.Ns * reg_x.s_, pk.DJpk2.Nss);
    DJ_Dist(v, pk.DJpk2, c_s_Ns_s);

    AddMod(v, v, PRF_ZZ(prf_key, pk.DJpk2.Ns), pk.DJpk2.Ns);

    if (b == 0){
        tau = - tau;
//        std::cout << "tau0: " << tau << std::endl;
        PowerMod(g_tau, pk.g, tau, pk.r_);
//        std::cout<<g_tau<<std::endl;
    } else {
//        g_tau = tau ;
//        std::cout << "tau1: " << tau << std::endl;
        PowerMod(g_tau, pk.g, tau, pk.r_);
//        std::cout<<g_tau<<std::endl;
    }
}

bool PVHSS_Ver(ZZ &v_ZZ, ZZ v0, ZZ v1, ZZ g_tau0, ZZ g_tau1, PVHSSPK pk, PVHSSPVK pvk) {

    ZZ right_side, left_side;

    v_ZZ = v1 - v0;

    PowerMod(right_side, pvk, v_ZZ, pk.r_);
    MulMod(left_side, g_tau0, g_tau1, pk.r_);

    if (left_side == right_side){
        std::cout << "Verification passed!\nThe value of y= " << v_ZZ << std::endl;
        return true;
    } else {
        printf("******************** ERROR ********************\n");
        std::cout << "Left side:" << left_side << std::endl;
        printf("\n\n");
        std::cout << "Right side:" << right_side << std::endl;
        std::cout << "The value of y= " << v_ZZ << std::endl;
        return false;
    }
}

