#include <time.h>
#include <sys/time.h>
#include "test.h"
// #include "Modular_Inversion.h"
#define num 1 // 0のとき平均クロックサイクル数を測定

int main(){
    mpz_t result,x,y,p,pre_com,new_pre_com,pre_hdfctmi,pre_BY;
	int loop,hloop,fltloop,byloop;
    mpz_init(result);
    mpz_init_set_str(x, "60", 10);
	mpz_init_set_str(y, "5", 10);
	mpz_t x1;
	mpz_init_set_str(x1, P_192_x1, 10);
	// mpz_and(x1,x1,x);
	// gmp_printf("x1=%Zd\n",x1);
	int check;
	check=384; // 192,224,256,384,512のいずれかを指定
	if(check==192){
		mpz_init_set_str(p, P_192, 10);
		mpz_init_set_str(pre_com, P_192_pre_com, 10);
		mpz_init_set_str(new_pre_com, P_192_new_pre_com, 10);
		loop = P_192_LOOPS_SICTMI;
		fltloop = P_192_LOOPS_FLT;
		byloop = P_192_LOOPS_BY;

	}else if(check==224){
		mpz_init_set_str(p, P_224, 10);
		mpz_init_set_str(pre_com, P_224_pre_com, 10);
		mpz_init_set_str(new_pre_com, P_224_new_pre_com, 10);
		mpz_init_set_str(pre_BY, P_224_PRE_BY, 10);
		loop = P_224_LOOPS_SICTMI;
		fltloop = P_224_LOOPS_FLT;
		byloop = P_224_LOOPS_BY;
	}else if(check==256){
		mpz_init_set_str(p, P_256, 10);
		mpz_init_set_str(pre_com, P_256_pre_com, 10);
		mpz_init_set_str(new_pre_com, P_256_new_pre_com, 10);
		loop = P_256_LOOPS_SICTMI;
		fltloop = P_256_LOOPS_FLT;
		byloop = P_256_LOOPS_BY;
	}else if(check==384){
		mpz_init_set_str(p, P_384, 10);
		mpz_init_set_str(pre_com, P_384_pre_com, 10);
		mpz_init_set_str(new_pre_com, P_384_new_pre_com, 10);
		mpz_init_set_str(pre_hdfctmi, P_384_PRE_hdFCIMI, 10);
		loop = P_384_LOOPS_SICTMI;
		hloop = P_384_LOOPS_hdFCTMI;
		fltloop = P_384_LOOPS_FLT;
		byloop = P_384_LOOPS_BY;
	}else if(check==512){
		mpz_init_set_str(p, P_512, 10);
		mpz_init_set_str(pre_com, P_512_pre_com, 10);
		mpz_init_set_str(new_pre_com, P_512_new_pre_com, 10);
		loop = P_512_LOOPS_SICTMI;
		fltloop = P_512_LOOPS_FLT;
		byloop = P_512_LOOPS_BY;
	}
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	gmp_randseed_ui(state, time(NULL));
	uint64_t CC1, CC2;


	//p以下、N個のランダムな整数を生成,test_num配列に格納
	int N = 100000;
	double per,a1,a2;
	mpz_t test_num[N];
	for(int i = 0; i < N; ++i){
		mpz_init(test_num[i]);
		mpz_urandomm(test_num[i], state, p);
	}

	if(num==0){
		// CC1 = GetCC();
		// for(int i = 0; i < N; ++i){
		// 	SICT_MI(result, test_num[i], p, pre_com, loop);
		// }
		// CC2 = GetCC();
		// 	printf("SICTMI\n%ld\n", (CC2-CC1)/N);
		// 	a1=(CC2-CC1)/N;
		// CC1 = GetCC();
		// for(int i = 0; i < N; ++i){
		// 	ex_SICT_MI(result,test_num[i], p, new_pre_com, loop-2);
		// }
		// CC2 = GetCC();
		// 	printf("卒論\n%ld\n", (CC2-CC1)/N);
		CC1 = GetCC();
		for(int i = 0; i < N; ++i){
			FLT(result,test_num[i], p, fltloop);
		}
		CC2 = GetCC();
			printf("FLT\n%ld\n", (CC2-CC1)/N);

		CC1 = GetCC();
		for(int i = 0; i < N; ++i){
			ECTMI2(result,test_num[i], p, pre_com, loop);
		}
		CC2 = GetCC();
			printf("JM23\n%ld\n", (CC2-CC1)/N);

		CC1 = GetCC();		
		for(int i = 0; i < N; ++i){
			KM1(result,test_num[i], p, new_pre_com, loop-2);
		}
		CC2 = GetCC();
			printf("KM1\n%ld\n", (CC2-CC1)/N);

		CC1 = GetCC();		
		for(int i = 0; i < N; ++i){
			KM2(&result,test_num[i], p, new_pre_com, loop-2);
		}
		CC2 = GetCC();
			printf("KM2\n%ld\n", (CC2-CC1)/N);
		CC1 = GetCC();		
		for(int i = 0; i < N; ++i){
			mpz_invert(result, test_num[i], p);
		}
		CC2 = GetCC();
			printf("KM2\n%ld\n", (CC2-CC1)/N);
		// per = 100*(1-a2/a1);
		// printf("割合:%f\n",per);
	}else{//逆元確認用
		// SICT_MI(result, x, p, pre_com, loop);
		// ex_SICT_MI(result, x, p, new_pre_com, loop-2);
		// FLT(result, x, p, fltloop);
		BY(result, x, p, pre_BY,byloop);
		// KM1(result, x, p, new_pre_com, loop-2);
		// ex_test2(result,x, p, new_pre_com, loop-2, x1);
		// ex_test3(&result, x, p, new_pre_com, loop-2);
		// ECTMI3(&result, x, p, new_pre_com, loop-2);
		mpz_invert(result, x, p);
		gmp_printf("逆元=%Zd\n",result);
	}

	return 0;

    mpz_clear(x); mpz_clear(p); mpz_clear(result);mpz_clear(pre_com);mpz_clear(y); 
    mpz_clear(new_pre_com);mpz_clear(x1);mpz_clear(pre_hdfctmi);mpz_clear(pre_BY);
}