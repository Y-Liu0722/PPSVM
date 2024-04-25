#pragma once
#include "helper.h"

typedef struct
{
    int s;
    ZZ N, Ns, Nss; // N, N^s, N^{s+1}
} DJ_PK;

void DJ_Gen(DJ_PK &pk, ZZ &d, int s);
void DJ_Enc(ZZ &ct, DJ_PK pk, ZZ x);
void DJ_Dec(ZZ &x, DJ_PK pk, ZZ sk, ZZ ct);

void DJ_Dist(ZZ &res, DJ_PK pk, ZZ x);
void DJ_Dist_Checker();

void DJ_exp(ZZ &res, DJ_PK pk, ZZ x); 
void DJ_log(ZZ &res, DJ_PK pk, ZZ x); 
bool DJ_Exp_Log_Checker();
