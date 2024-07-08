#include "PKE.h"

void SetPara(PKE_Para & pkePara) {
    switch(pkePara.msg_bit){
        case 2:
            pkePara.N=8192;
            pkePara.p_bit=63;
            pkePara.q_bit=148;
            break;
        case 16:
            pkePara.N=8192;
            pkePara.p_bit=77;
            pkePara.q_bit=176;
            break;
        case 32:
            pkePara.N=8192;
            pkePara.p_bit=93;
            pkePara.q_bit=208;
            break;
        case 64:
            pkePara.N=16384;
            pkePara.p_bit=126;
            pkePara.q_bit=275;
            break;
        case 128:
            pkePara.N=16384;
            pkePara.p_bit=190;
            pkePara.q_bit=403;
            break;
        case 256:
            pkePara.N=32768;
            pkePara.p_bit=319;
            pkePara.q_bit=662;
            break;
        case 512:
            pkePara.N=65536;
            pkePara.p_bit=576;
            pkePara.q_bit=1177;
            break;
        case 1024:
            pkePara.N=131072;
            pkePara.p_bit=1088;
            pkePara.q_bit=2204;
            break;
        default:printf("error\n"); break;
    }
}

void PKE_Gen(PKE_Para &pkePara, PKE_PK &pkePk, PKE_SK &pkeSk) {
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    //initialize the parameters
    SetPara(pkePara);

    RandomPrime(pkePara.r, pkePara.msg_bit);
    while (pkePara.r % (2 * pkePara.N) != 1){
        RandomPrime(pkePara.r, pkePara.msg_bit);
    }

//    std::cout << "r: " << pkePara.r << std::endl;

    NTL::power(pkePara.p, 2, pkePara.p_bit);
    NTL::power(pkePara.q, 2, pkePara.q_bit);
//    NTL::mul(pkePara.q, pkePara.q, pkePara.p);
    pkePara.r_context = ZZ_pContext(pkePara.r);
    pkePara.q_context = ZZ_pContext(pkePara.q);
    ZZ_p::init(pkePara.q);
    SetCoeff(pkePara.xN, 0, 1);
    SetCoeff(pkePara.xN, pkePara.N, 1);

    {
        // Switch to modulus r
        ZZ_pPush push(pkePara.r_context);
        build(pkePara.modulus_r, pkePara.xN);
//        pkePara.num_slots = 1024;
        RootEDF(pkePara.factors, pkePara.modulus_r);

//        pkePara.d = pkePara.N / pkePara.num_slots;

        vec_ZZ_pX n_star;
        vec_ZZ_pX t;

        pkePara.modulus_f.SetLength(pkePara.N);
        n_star.SetLength(pkePara.N);
        t.SetLength(pkePara.N);
        pkePara.e.SetLength(pkePara.N);
        for (auto i = 0; i < pkePara.N; i++){
            build(pkePara.modulus_f[i], pkePara.factors[i]);
            n_star[i] = pkePara.xN / pkePara.factors[i];
            InvMod(t[i], n_star[i] % pkePara.factors[i], pkePara.factors[i]);
            pkePara.e[i] = t[i] * n_star[i];
        }
    }
//    std::cout << "good" << std::endl;

    pkePara.twice_p = 2 * pkePara.p;
    pkePara.twice_q = 2 * pkePara.q;
    pkePara.half_p = pkePara.p / 2;


//    std::cout << "gen sk" << std::endl;
//    ZZ_pXModulus modulus(pkePara.xN);
    build(pkePara.modulus, pkePara.xN);
    //gen sk
    ZZ_pX hat_s, e;
    Random_ZZ_pX(pkePk[0], pkePara.N, pkePara.q_bit);
    SecretKey(hat_s, pkePara.N, pkePara.hsk);
    GaussRand(e, pkePara.N);

    MulMod(pkePk[1], pkePk[0], hat_s, pkePara.modulus);
    pkePk[1] = pkePk[1] + e;
    //gen sk
    SetCoeff(pkeSk[0], 0, 1);
    pkeSk[1] = hat_s;

    conv(pkePara.q_div_p, pkePara.q / pkePara.p);
}

void PKE_Enc(PKE_Ciphertext &c, const PKE_Para& pkePara, const PKE_PK &pkePk, const ZZ_pX &m) {
    ZZ_pX v, e1, e2;
//    ZZ temp = pkePara.q / pkePara.p;
//    ZZ_p q_div_p;
//    conv(q_div_p, temp);
    SecretKey(v, pkePara.N, pkePara.hsk);
    GaussRand(e1, pkePara.N);
    GaussRand(e2, pkePara.N);

    MulMod(c[0], pkePk[1], v, pkePara.modulus);

    c[0] = c[0] + e2 + pkePara.q_div_p * m;
    MulMod(c[1], pkePk[0], v, pkePara.modulus);
    c[1] = e1 - c[1];
}

void PKE_OKDM(Ciphertext &C, const PKE_Para &pkePara, const PKE_PK &pkePk, const ZZ_pX &m) {
//    PKE_Ciphertext c_m, c_ms;
    C.c_m.SetLength(2);
    C.c_ms.SetLength(2);

    PKE_Enc(C.c_m, pkePara, pkePk, m);
    PKE_Enc(C.c_ms, pkePara, pkePk, ZZ_pX::zero());

    C.c_ms[1] = C.c_ms[1] + pkePara.q_div_p * m;
}

void PKE_DDec(MemoryV &t_xy, const PKE_Para& pkePara, MemoryV t_x, Ciphertext C_y) {
    t_xy.SetLength(2);
    ZZ coeff;
    ZZX temp;
    ZZ_pX temp1, temp2;
    MulMod(temp1, t_x[0], C_y.c_m[0], pkePara.modulus);
    MulMod(temp2, t_x[1], C_y.c_m[1], pkePara.modulus);
    temp1 = temp1 + temp2;
    conv(temp,temp1);
    for (int i = 0; i < pkePara.N; i++) {
        GetCoeff(coeff, temp, i);
        coeff = (coeff * pkePara.twice_p + pkePara.q) / (pkePara.twice_q);
        coeff= coeff % pkePara.p;
        if(coeff > pkePara.half_p){
            coeff-=pkePara.p;
        }
        SetCoeff(temp, i, coeff);
    }
    conv(t_xy[0],temp);
    MulMod(temp1, t_x[0], C_y.c_ms[0], pkePara.modulus);
    MulMod(temp2, t_x[1], C_y.c_ms[1], pkePara.modulus);
    temp1 = temp1 + temp2;
    conv(temp,temp1);
    for (int i = 0; i < pkePara.N; i++) {
        GetCoeff(coeff, temp, i);
        coeff = (coeff * pkePara.twice_p + pkePara.q) / (pkePara.twice_q);
        coeff= coeff % pkePara.p;
        if(coeff > pkePara.half_p){
            coeff-=pkePara.p;
        }
        SetCoeff(temp, i, coeff);
    }
    conv(t_xy[1],temp);
}

// CRT encoding
ZZ_pX Encode(const PKE_Para& pkePara, const vec_ZZ_p &x){
    ZZ_pPush push(pkePara.r_context);
//    clear(hat_x);
    vec_ZZ_p x_(x);
    ZZ_pX hat_x;
    if (x_.length() > pkePara.N){
        std::cerr << "Encode not support vector with length > num of slots" << std::endl;
    }
    if (x.length() < pkePara.N){
        x_.SetLength(pkePara.N, ZZ_p::zero());
    }
//    ZZ_pX temp;
    for (auto i = 0; i < pkePara.N; i++){
//        SetCoeff(temp, 0, x_[i]);
        hat_x = hat_x + x_[i] * pkePara.e[i];
    }
    hat_x = hat_x % pkePara.modulus_r;
//    ZZ_p::init(pkePara.q);
    return hat_x;
}

// CRT decoding
vec_ZZ_p Decode(const PKE_Para& pkePara, const ZZ_pX& hat_x){
//    ZZ_p::init(pkePara.r);
    vec_ZZ_p x;
    x.SetLength(pkePara.N);
//    ZZ_pContext context;
//    context.save();
//    NTL_EXEC_RANGE(pkePara.N, first, last)
////        context.restore();
//        ZZ_pPush push(pkePara.r_context);
        ZZ_pX px;
        // px.SetLength(pkePara.N);
//        std::cout << first << std::endl;
//        std::cout << last << std::endl;
        for (auto i = 0; i < pkePara.N; i++) {
            // std::cout << pkePara.factors[i]
            rem(px, hat_x, pkePara.modulus_f[i]);
            GetCoeff(x[i], px, 0);
        }
//    NTL_EXEC_RANGE_END
//    ZZ_p::init(pkePara.q);
    return x;
//    std::cout << std::endl;
}
