#include <time.h>
#include <sys/time.h>
#include "test.h"
// #include "Modular_Inversion.h"
#define num 0 // 0のとき平均クロックサイクル数を測定


int main(){

	mpz_t temp1; mpz_init(temp1);
	mpz_t temp2; mpz_init(temp2);
	mpz_t temp3; mpz_init(temp3);
	mpz_t temp4; mpz_init(temp4);
	mpz_t *sort1[2]; sort1[0] = &temp1; sort1[1] = &temp2;
	mpz_t *sort2[2]; sort2[0] = &temp3; sort2[1] = &temp4;
	bool b0,b1;
	b0=0;
	b1=1;

    mpz_t result,x,y,p,pre_com,new_pre_com,pre_hdfctmi;
	int loop,hloop,fltloop;
    mpz_init(result);
    mpz_init_set_str(x, "0", 10);
	mpz_init_set_str(y, "1", 10);
	mpz_t x1;
	mpz_init_set_str(x1, P_192_x1, 10);
	int check;
	check=192; // 192,224,256,384,512のいずれかを指定
	if(check==192){
		mpz_init_set_str(p, P_192, 10);
		mpz_init_set_str(pre_com, P_192_pre_com, 10);
		mpz_init_set_str(new_pre_com, P_192_new_pre_com, 10);
		loop = P_192_LOOPS_SICTMI;
		fltloop = P_192_LOOPS_FLT;
	}else if(check==224){
		mpz_init_set_str(p, P_224, 10);
		mpz_init_set_str(pre_com, P_224_pre_com, 10);
		mpz_init_set_str(new_pre_com, P_224_new_pre_com, 10);
		loop = P_224_LOOPS_SICTMI;
		fltloop = P_224_LOOPS_FLT;
	}else if(check==256){
		mpz_init_set_str(p, P_256, 10);
		mpz_init_set_str(pre_com, P_256_pre_com, 10);
		mpz_init_set_str(new_pre_com, P_256_new_pre_com, 10);
		loop = P_256_LOOPS_SICTMI;
		fltloop = P_256_LOOPS_FLT;
	}else if(check==384){
		mpz_init_set_str(p, P_384, 10);
		mpz_init_set_str(pre_com, P_384_pre_com, 10);
		mpz_init_set_str(new_pre_com, P_384_new_pre_com, 10);
		loop = P_384_LOOPS_SICTMI;
		fltloop = P_384_LOOPS_FLT;
	}else if(check==512){
		mpz_init_set_str(p, P_512, 10);
		mpz_init_set_str(pre_com, P_512_pre_com, 10);
		mpz_init_set_str(new_pre_com, P_512_new_pre_com, 10);
		loop = P_512_LOOPS_SICTMI;
		fltloop = P_512_LOOPS_FLT;
	}
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	gmp_randseed_ui(state, time(NULL));
	uint64_t CC1, CC2;
	uint64_t Mp, Sp, Ap, Subp, rshift, lshift, select,and_p,or_p, xor_p, cmp, bitflip;
	uint64_t Mpw, andpw, orpw, Apw, orw;
	uint64_t andw, xorw, Subw, rshiftw, lshiftw;


	//p以下、N個のランダムな整数を生成,test_num配列に格納
	int N = 200000;
	double per,a1,a2;
	mpz_t test_num[N];
	mpz_t test_num2[N];
	// mpz_t test_num3[N];
	// mpz_t W;
	// mpz_init(W);
	for(int i = 0; i < N; ++i){
		mpz_init(test_num[i]);
		mpz_urandomm(test_num[i], state, p);
		mpz_init(test_num2[i]);
		mpz_urandomm(test_num2[i], state, p);
	}

	// for(int i = 0; i < N; ++i){
	// 	mpz_init(test_num[i]);
	// 	mpz_urandomm(test_num3[i], state, p);
	// }

	// for(int i = 0; i < 5; ++i){
	// 	gmp_printf("%Zd\n",test_num[i]);
	// }
	// printf("\n");
	// for(int i = 0; i < 5; ++i){
	// 	gmp_printf("%Zd\n",test_num2[i]);
	// }

	if(num==0){
		CC1 = GetCC();
		for(int i = 0; i < N; ++i){
			mpz_mul(result, test_num[i], test_num2[i]);
			mpz_mod(result, result, p);
		}
		CC2 = GetCC();
		Mp= (CC2-CC1)/N;
		CC1 = GetCC();
		for(int i = 0; i < N; ++i){
			mpz_mul(result, test_num[i], test_num[i]);
			mpz_mod(result, result, p);
		}
		CC2 = GetCC();
		Sp= (CC2-CC1)/N;
		CC1 = GetCC();
		for(int i = 0; i < N; ++i){
			mpz_add(result, test_num[i], test_num2[i]);
			// mpz_add(result, test_num2[i], test_num[i]);
		}
		CC2 = GetCC();
		Ap = (CC2-CC1)/N;
		CC1 = GetCC();		
		for(int i = 0; i < N; ++i){
			mpz_sub(result, test_num[i], test_num2[i]);
		}
		CC2 = GetCC();
		Subp = (CC2-CC1)/N;
		CC1 = GetCC();		
		for(int i = 0; i < N; ++i){
			mpz_mul_2exp(result, test_num[i], 1);
		}
		CC2 = GetCC();
		rshift = (CC2-CC1)/N;
		CC1 = GetCC();		
		for(int i = 0; i < N; ++i){
			mpz_tdiv_q_2exp(result, test_num[i], 1);
		}
		CC2 = GetCC();
		lshift = (CC2-CC1)/N;
		CC1 = GetCC();		
		for(int i = 0; i < N; ++i){
			mpz_set(test_num[i], *sort1[b1]);
			mpz_set(test_num2[i], *sort1[!b1]);
			mpz_set(test_num2[i], *sort2[b1]); 
			mpz_set(test_num[i], *sort2[!b1]);
		}
		CC2 = GetCC();
		select = (CC2-CC1)/N;
		CC1 = GetCC();		
		for(int i = 0; i < N; ++i){
			mpz_xor(result, test_num[i], y);
		}
		CC2 = GetCC();
		xor_p = (CC2-CC1)/N;
		CC1 = GetCC();		
		for(int i = 0; i < N; ++i){
			mpz_cmp(test_num[i], test_num2[i]);
		}
		CC2 = GetCC();
		cmp = (CC2-CC1)/N;
		CC1 = GetCC();		
		for(int i = 0; i < N; ++i){
			mpz_com(result, test_num[i]);
		}
		CC2 = GetCC();
		bitflip = (CC2-CC1)/N;

		CC1 = GetCC();
		for(int i = 0; i < N; ++i){
			mpz_mul_ui(result, test_num[i], b0);
			mpz_mul_ui(result, test_num[i], b1);
		}
		CC2 = GetCC();
		Mpw= (CC2-CC1)/(2*N);
		CC1 = GetCC();
		for(int i = 0; i < N; ++i){
			mpz_and(result, test_num[i], x);
		}
		CC2 = GetCC();
		andpw= (CC2-CC1)/N;
		CC1 = GetCC();
		for(int i = 0; i < N; ++i){
			mpz_ior(result, test_num[i], x);
		}
		CC2 = GetCC();
		orpw= (CC2-CC1)/N;

		CC1 = GetCC();
		for(int i = 0; i < N; ++i){
			mpz_and(result, test_num[i], test_num2[i]);
		}
		CC2 = GetCC();
		and_p= (CC2-CC1)/N;
		CC1 = GetCC();
		for(int i = 0; i < N; ++i){
			mpz_ior(result, test_num[i], test_num2[i]);
		}
		CC2 = GetCC();
		or_p= (CC2-CC1)/N;
		CC1 = GetCC();
		for(int i = 0; i < N; ++i){
			mpz_add(result, test_num[i], y);
		}
		CC2 = GetCC();
		Apw= (CC2-CC1)/N;
		CC1 = GetCC();
		for(int i = 0; i < N; ++i){
			b0 || b1;
		}
		CC2 = GetCC();
		orw= (CC2-CC1)/N;

		CC1 = GetCC();
		for(int i = 0; i < N; ++i){
			b0 & b1;
		}
		CC2 = GetCC();
		andw= (CC2-CC1)/N;

		CC1 = GetCC();
		for(int i = 0; i < N; ++i){
			b0 & b1;
		}
		CC2 = GetCC();
		andw= (CC2-CC1)/N;

		CC1 = GetCC();
		for(int i = 0; i < N; ++i){
			b0 ^ b1;
		}
		CC2 = GetCC();
		xorw= (CC2-CC1)/N;

		CC1 = GetCC();
		for(int i = 0; i < N; ++i){
			b1-b0;
		}
		CC2 = GetCC();
		Subw= (CC2-CC1)/N;

		CC1 = GetCC();
		for(int i = 0; i < N; ++i){
			b1<<1;
			b0<<1;
		}

		CC2 = GetCC();
		rshiftw= (CC2-CC1)/(2*N);

		printf("Mp,Sp,Ap,Subp,select,rshift,lshift,xorp,andp,bitflipp,cmp,Mpw,andpw,orpw\n");
		printf("%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld\n", Mp, Sp, Ap, Subp, select, rshift,lshift,xor_p, and_p, bitflip,cmp, Mpw, andpw, orpw);
		// // printf("Sp/Mp,Ap/Mp,Subp/Mp,rshift/Mp,lshift/Mp,xorp/Mp,and_p/Mp,bitflip/Mp,cmp/Mp,Mpw/Mp,andpw/Mp,orpw/Mp\n");
		// printf("Sp,Ap,Subp,select,rshift,lshift,xorp,andp,bitflipp,cmp,Mpw,andpw,orpw\n");
		// printf("%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n",   (double)Sp/Mp, (double)Ap/Mp, (double)Subp/Mp, (double)select/Mp, (double)rshift/Mp, (double)lshift/Mp, (double)xor_p/Mp, (double)and_p/Mp, (double)bitflip/Mp, (double)cmp/Mp, (double)Mpw/Mp, (double)andpw/Mp, (double)orpw/Mp);
		printf("Apw, orw\n");
		printf("%ld,%ld\n",Apw, orw);
		// printf("%.2f,%.2f\n", (double)Apw/Mp,(double)orw/Mp);

		printf("andw, xorw, Subw, rshiftw\n");
		printf("%ld,%ld,%ld,%ld\n",andw, xorw, Subw, rshiftw);
		// printf("%.2f,%.2f,%.2f\n", (double)andw/Mp,(double)xorw/Mp,(double)Subw/Mp);

		// printf("Mpw/Mp, and_p/Mp, or_p/Mp, andpw/Mp, orpw/Mp\n");
		// printf("%f,%f,%f,%f,%f\n",(double)Mpw/Mp,(double)and_p/Mp, (double)or_p/Mp, (double)andpw/Mp, (double)orpw/Mp);
		// printf("andw\n%ld\n", (CC2-CC1)/N);
	}else{//逆元確認用
		// SICT_MI(result, x, p, pre_com, loop);
		// ex_SICT_MI(result, x, p, new_pre_com, loop-2);
		// KM1(result, x, p, new_pre_com, loop-2);
		// ex_test2(result,x, p, new_pre_com, loop-2, x1);
		// ex_test3(&result, x, p, new_pre_com, loop-2);
		// ECTMI3(&result, x, p, new_pre_com, loop-2);
		// mpz_invert(result, x, p);
		// gmp_printf("逆元=%Zd\n",result);
	}

	return 0;

    mpz_clear(x); mpz_clear(p); mpz_clear(result);mpz_clear(pre_com);mpz_clear(y); 
    mpz_clear(new_pre_com);mpz_clear(x1);
	mpz_clear(temp1); mpz_clear(temp2);
	mpz_clear(temp3); mpz_clear(temp4);
	
}