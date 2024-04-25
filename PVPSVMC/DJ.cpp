#include "DJ.h"

#ifndef INVI
#define INVI

Vec<ZZ> inv_i;
Vec<ZZ> exp_N_i;
Vec<ZZ> log__N_i_1;

#endif

void DJ_Gen(DJ_PK &pk, ZZ &sk, int s) {
    ZZ p, q;
    pk.s = s;
    // 128 bit
//        GenGermainPrime(p, 1536); // safe prime
//        GenGermainPrime(q, 1536);

    p = conv<ZZ>(
            "241031242692103258858011660602831411291209324794568895135967503906525739"
            "159180320066908502410734604966344876628088800478786241697879495832496961"
            "298789077465145521333938162522477078207791768149967684554313738782005759"
            "734585790459910946138712209950796499781564134230067762947335528161742841"
            "179416396778587037036896910922159194305423201156275845008057958785090099"
            "371489228347664663118151506380487337518226050624699283789870597101252584"
            "3324401232986857004760339321639");
    q = conv<ZZ>(
            "241031242692103258858011660602831411291209324794568895135967503906525739"
            "159180320066908502410734604966344876628088800478786241697879495832496961"
            "298789077465145521333938162522477078207791768149967684554313738782005759"
            "734585790459910946138712209950796499781564134230067762947335528161742841"
            "179416396778587037036896910922159194305423201156275845008057958785090099"
            "371489228347664663118151506380487337518226050624699283789870597101252584"
            "3324401232986857004760339319223");

    // 118 bit
    // GenGermainPrime(p, 1024); // safe prime
    // GenGermainPrime(q, 1024);
    // p =
    // conv<ZZ>("166233368868847378150709182712632437464545363054158394531844328754558031827000690649980479925587807047206364466579401042695955597766505153469629285725614151722406341121412138060966963852833451734571031793524529752601181071118077244607107664780951337127495029527362834055643232430411000382851603611712401019243");
    // q =
    // conv<ZZ>("101880783786260223305427009662645830161677199493002793826918563945056325203550598774895893664897969187660499167834545419951427565438111153189836460715802235978731557523542988036414173096593856568950867698439044589562662387659086154436677345217663111722570654601519389883320023346908316033140215862510577762751");

    mul(pk.N, p, q);
    power(pk.Ns, pk.N, pk.s);
    mul(pk.Nss, pk.Ns, pk.N);
    LCM(sk, p - 1, q - 1); // phi(N)

    ZZ inv = conv<ZZ>(1);
    ZZ minus = conv<ZZ>(-1);
    for (auto i = 2; i <= s; i++) {
        ZZ inv_i;
        InvMod(inv_i, ZZ(i), pk.Nss);
        mul(inv, inv, inv_i);
        ZZ tmp;
        PowerMod(tmp, pk.N, i, pk.Nss);
        MulMod(tmp, tmp, inv, pk.Nss);
        exp_N_i.append(tmp);

        ZZ minus_i;
        power(minus_i, minus, i - 1);
        PowerMod(tmp, pk.N, i - 1, pk.Ns);
        mul(tmp, tmp, minus_i);
        MulMod(tmp, tmp, inv_i, pk.Ns);
        log__N_i_1.append(tmp);
    }

//    inv_i.SetLength(2);
//    InvMod(inv_i[0], ZZ(2), pk.Nss);
//    InvMod(inv_i[1], ZZ(3), pk.Nss);
}

void DJ_Enc(ZZ &ct, DJ_PK pk, ZZ x) {
    ZZ r, temp1, temp2;

    RandomBnd(r, pk.Nss);
    PowerMod(temp1, r, pk.Ns, pk.Nss);
    DJ_exp(temp2, pk, x);
    MulMod(ct, temp1, temp2, pk.Nss);
}

void DJ_Dec(ZZ &x, DJ_PK pk, ZZ sk, ZZ ct) {
    ZZ temp, inv_sk;
    PowerMod(temp, ct, sk, pk.Nss);
    DJ_log(temp, pk, temp);
    InvMod(inv_sk, sk, pk.Nss);

    MulMod(x, temp, inv_sk, pk.Ns);
}

void DJ_exp(ZZ &res, DJ_PK pk, ZZ x) {
    ZZ Nx;
    MulMod(Nx, pk.N, x, pk.Nss);
    AddMod(res, ZZ(1), Nx, pk.Nss); // 1 + Nx

    ZZ tmp = Nx;
    for (int i = 2; i <= pk.s; ++i) {
        power(tmp, x, i);
        MulMod(tmp, tmp, exp_N_i[i-2], pk.Nss);
//        MulMod(tmp, tmp, Nx, pk.Nss);
//        MulMod(tmp, tmp, inv_i[i - 2], pk.Nss);
        AddMod(res, res, tmp, pk.Nss);
    }
}

void DJ_log(ZZ &res, DJ_PK pk, ZZ x) {
    ZZ newx;
    div(newx, x - 1, pk.N);
    res = newx;
    ZZ tmp;
//    ZZ tmp = newx;
//    ZZ minusNx;
//    MulMod(minusNx, pk.N, -newx, pk.Ns);
    for (int i = 2; i <= pk.s; ++i) {
        power(tmp, newx, i);
        MulMod(tmp, tmp, log__N_i_1[i-2], pk.Ns);

//        MulMod(tmp, tmp, minusNx, pk.Ns);
//        MulMod(tmp, tmp, inv_i[i - 2], pk.Ns);
        AddMod(res, res, tmp, pk.Ns);
//        MulMod(tmp, tmp, ZZ(i), pk.Ns);
    }
}

bool DJ_Exp_Log_Checker() {
    DJ_PK pk;
    ZZ sk;
    DJ_Gen(pk, sk, 3);

    ZZ m;
    ZZ ct, dct;

    RandomBnd(m, pk.Ns);

    ZZ n;
    RandomBnd(n, pk.N);
    ZZ expm, expn, expmn1, expmn2;
    DJ_exp(expm, pk, m);
    DJ_exp(expn, pk, n);
    MulMod(expmn1, expm, expn, pk.Nss);
    DJ_exp(expmn2, pk, m + n);
    if (expmn1 != expmn2) {
        return false;
    }

    ZZ logm;
    DJ_log(logm, pk, m * pk.N + 1);
    DJ_exp(expm, pk, logm);
    if (pk.N * m + 1 != expm) {
        return false;
    }

    DJ_exp(expm, pk, m);
    DJ_log(logm, pk, expm);
    if (logm != m) {
        return false;
    }

    DJ_Enc(ct, pk, m);
    DJ_Dec(dct, pk, sk, ct);

    if (m != dct) {
        return false;
    }
    return true;
}

void DJ_Dist(ZZ &res, DJ_PK pk, ZZ c) {
    // c mod N, (c mod N)^-1, \frac{c}{c mod N}
    ZZ c_N, inv_c_N, c_c_N;
    rem(c_N, c, pk.N);

//    std::cout << "c: "<< c << std::endl;
//    std::cout << "c mod N: " << c_N << std::endl;
//    std::cout << "c / (c mod N): " << c/c_N << std::endl;
//    std::cout << "c / (c mod N) *  (c mod N): " << (c/c_N) * c_N << std::endl;
//    std::cout << "diff: " << (c/c_N) * c_N - c << std::endl;

    ZZ res_c, res_c_N;
//    DJ_log(res_c, pk, c);
//    DJ_log(res_c_N, pk, c_N);
    //    AddMod(c_N, c, ZZ(0), pk.N);

    InvMod(inv_c_N, c_N, pk.Nss);


    MulMod(c_c_N, c, inv_c_N, pk.Nss);

    DJ_log(res, pk, c_c_N);
}

void DJ_Dist_Checker() {
    DJ_PK pk;
    ZZ sk;
    DJ_Gen(pk, sk, 2);

    ZZ c, x;
    RandomBnd(c, pk.Nss);
    RandomBnd(x, pk.Ns);

    // c exp(x)
    ZZ ex, cex;
    DJ_exp(ex, pk, x);
    MulMod(cex, c, ex, pk.Nss);

    //    ZZ tmp;
    //    rem(tmp,  ex, pk.N);
    //    std::cout << tmp << std::endl;
    //
    //    rem(tmp,  c, pk.N);
    //    std::cout << tmp << std::endl;
    //
    //    rem(tmp, cex, pk.N);
    //    std::cout << tmp << std::endl;

    ZZ res1, res2;
    DJ_Dist(res1, pk, cex);
    DJ_Dist(res2, pk, c);
    std::cout << res1 << std::endl;
    std::cout << res2 << std::endl;

    ZZ x_;
    SubMod(x_, res1, res2, pk.Ns);

    std::cout << x_ << std::endl;
    std::cout << x << std::endl;
    std::cout << (x == x_) << std::endl;
}
