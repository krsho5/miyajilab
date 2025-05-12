#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <gmp.h>

void pornin(mpz_t m, mpz_t y, mpz_t k, mpz_t len){
mpz_t a,u,b,v;
mpz_init_set(a, y);mpz_init_set(b, m);mpz_init_set_str(u, "1", 10);mpz_init_set_str(v, "0", 10);
mpz_t l,tempk;
mpz_init(l);
mpz_sub(tempk,k,u);//tempk=k-1
mpz_cdiv_q(l,len, tempk);//l:ループ回数
long loop;
loop = mpz_get_si(l);//mpz_tからunsigned long
for(long i = 0; i < loop; ++i){
    size_t len_a,len_b,K;//K=2k
    len_a = mpz_sizeinbase(a,2);
    len_b = mpz_sizeinbase(b,2);
    K = mpz_get_si(k);
    // mpz_t len_A,len_B;
    // mpz_set_ui(len_A, (unsigned long)len_a);
    // mpz_set_ui(len_B, (unsigned long)len_b);
    long n;
    n = len_b;
    if(len_a > n){
        n = len_a;
    }
    if(K > n){
        n = K; 
    }
    mpz_t false_a,false_b;
    


    
}
mpz_clear(a);mpz_clear(u);mpz_clear(b);mpz_clear(v);
}