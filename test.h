#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <gmp.h>

/*SICTMI: The number of loops defined by 2*bitlen(p) with p224, p256, p384.*/
#define P_192_LOOPS_SICTMI 384
#define P_224_LOOPS_SICTMI 448 
#define P_256_LOOPS_SICTMI 512
#define P_384_LOOPS_SICTMI 768
#define P_512_LOOPS_SICTMI 1024

/*The number of pre_com 2^(-l) mod p*/
#define P_192_pre_com "680564733841876926908302470789826871294"
#define P_224_pre_com  "18831305204698540654176465351577770740545525012632715657217"
#define P_256_pre_com  "161759680002903838760694582340907787320905632487091026015073739997175"
#define P_384_pre_com  "39401997719623594004609818748773390158923402225834345130666249913583939917682844212927918355281735000658345862363656"
#define P_512_pre_com "104174602983677100090188404798997345435175366322587289367637992126922563506503570441557461874283519919404476999295620665941816345906030878912480330037339671"

/*The number of pre_com 2^(-l-2) mod p*/
#define P_192_new_pre_com "2722258935367507707633209883159307485176"
#define P_224_new_pre_com "75325220818794162616705861406311082962182100050530862628868"
#define P_256_new_pre_com "647038720011615355042778329363631149283622529948364104060294959988700"
#define P_384_new_pre_com "39401972289310938381602154694662719220454391091941040518820119441598594356240765864569875156350125997048562530117667"
#define P_512_new_pre_com "14464174036430487373532869249813997916320490672577356138845125196037333123807873035824224602633369394547812422623024511206003735188283805810620373131106399"
/*Definition of the field.*/
#define P_192 "6277101735386680763835789423207666416083908700390324961279"
#define P_224 "26959946667150639794667015087019630673557916260026308143510066298881"
#define P_256 "115792089210356248762697446949407573530086143415290314195533631308867097853951"
#define P_384 "39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319"
#define P_512 "134078079299425970995740249982058461274793658205923933777235614437217640300735469576801874298166903427690031858186486050853753882811946569946433649006084095"

/*Precomputed numbers of hdFCTMI.*/
#define P_224_PRE_hdFCIMI "3369993333197670545056871185872896170421996906849403135584645414912"
#define P_256_PRE_hdFCIMI "49471717823667788566512707997982133781718357895554504641135267859606011904"
#define P_384_PRE_hdFCIMI "36092402410143488406830006651847393538951028178043571563308995486226118841864038653366404110975726264988796028928"

/*hdFCTMI: The number of loops when d = 224, 256, or 384.*/
#define P_224_LOOPS_hdFCTMI 517 
#define P_256_LOOPS_hdFCTMI 590
#define P_384_LOOPS_hdFCTMI 885

/*2^(l+1) - 1*/
#define P_192_x1 "9850501549098619803069760025035903451269934817616361666987073351061430442874302652853566563721228910201656997576703"
void SICT_MI(mpz_t v, mpz_t a, mpz_t p, mpz_t pre_com,int LOOPS){
	mpz_t u,q,r;
	mpz_init_set(u, a); mpz_set(v, p);
	mpz_init_set_str(q, "0", 10);
	mpz_init_set_str(r, "1", 10);

	mpz_t temp1; mpz_init(temp1);
	mpz_t temp2; mpz_init(temp2);
	mpz_t temp3; mpz_init(temp3);
	mpz_t temp4; mpz_init(temp4);
	bool u_lsb, v_lsb, b1;
	
	mpz_t *sort1[2]; sort1[0] = &temp1; sort1[1] = &temp2;
	mpz_t *sort2[2]; sort2[0] = &temp3; sort2[1] = &temp4;
	int result;

	for (int i = 0; i < LOOPS; ++i){
		mpz_sub(temp1, v, u); //sub v-u
		mpz_sub(temp3, q, r); //sub q-r
		u_lsb = mpz_tstbit(u, 0);
		v_lsb = mpz_tstbit(v, 0);
		
		b1 = (!u_lsb) && v_lsb;
		//mpz_set(temp2, u); temp2->_mp_size *= b1;
		mpz_mul_ui(temp2, u, b1);
		//mpz_set(temp4, r); temp4->_mp_size *= b1;
		mpz_mul_ui(temp4, r, b1);

		b1 = u_lsb && (!v_lsb);
		//v->_mp_size *= b1;
		mpz_mul_ui(v, v, b1);
		mpz_xor(temp2, temp2, v);//xor
		//q->_mp_size *= b1;
		mpz_mul_ui(q, q, b1);
		mpz_xor(temp4, temp4, q);//xor

		b1 = u_lsb && v_lsb;
		//mpz_set(v, temp1); v->_mp_size *= b1;
		mpz_mul_ui(v, temp1, b1);
		mpz_xor(temp2, temp2, v);//xor
		mpz_tdiv_q_2exp(temp2, temp2, 1);//rightshift
		//u->_mp_size *= b1;
		mpz_mul_ui(u, u, b1);
		//mpz_set(q, temp3); q->_mp_size *= b1;
		mpz_mul_ui(q, temp3, b1);
		mpz_xor(temp4, temp4, q);//xor
		//r->_mp_size *= b1;
		mpz_mul_ui(r, r, b1);

		b1 = (!u_lsb) || (!v_lsb);
		//temp1->_mp_size *= b1;
		mpz_mul_ui(temp1, temp1, b1);
		mpz_xor(temp1, temp1, u);//xor
		//temp3->_mp_size *= b1;
		mpz_mul_ui(temp3, temp3, b1);
		mpz_xor(temp3, temp3, r);//xor
		mpz_mul_2exp(temp3, temp3, 1);//leftshift

        //gmp_printf("(temp1,temp2,temp3,temp4)=(%Zd,%Zd,%Zd,%Zd)\n",temp1,temp2,temp3,temp4);
		
		
		result = mpz_cmp(temp2, temp1);
		if (result >= 0) {
			b1 = 1;
			//printf("z=%d\n",result);
		} else {
			b1 = 0;
			//printf("z=%d\n",result);
		}

		mpz_set(v, *sort1[b1]); mpz_set(u, *sort1[!b1]);
		mpz_set(q, *sort2[b1]); mpz_set(r, *sort2[!b1]);
		// printf("%d  ",i);
		// gmp_printf("(u,v,q,r)=(%Zd,%Zd,%Zd,%Zd)\n",u,v,q,r);
		// if(mpz_cmp(u,x)==0){
		// 	// gmp_printf("%Zd\n",u);
		// 	printf("%d,",i);
		// 	break;
		// }
	}
	mpz_mul(q, q, pre_com);	mpz_mod(q, q, p);
	// printf("%d\n",sum/LOOPS);
	// gmp_printf("逆元=%Zd\n", q);
	mpz_clear(u); mpz_clear(temp1); mpz_clear(temp2); 
	mpz_clear(q); mpz_clear(r); mpz_clear(temp3); mpz_clear(temp4); 
	return ;
}


// 0<=y<m m奇数
void ex_SICT_MI(mpz_t v, mpz_t y, mpz_t m, mpz_t pre_com,int LOOPS){
	mpz_t a,b,u;
	mpz_init_set(a, y);
	mpz_init_set(b, m);
	mpz_init_set_str(u, "1", 10);
	mpz_init_set_str(v, "0", 10);
	mpz_t tempa; mpz_init(tempa);
	mpz_t tempb; mpz_init(tempb);
	mpz_t tempu; mpz_init(tempu);
	mpz_t tempv; mpz_init(tempv);
	mpz_t tempa1; mpz_init(tempa1);
	mpz_t tempu1; mpz_init(tempu1);
	bool s, z, b1;
	int result;
	//gmp_printf("(a,u,b,v)=(%Zd,%Zd,%Zd,%Zd)\n",a,u,b,v);
	for (int i = 0; i < LOOPS; ++i){
		// printf("%d\n",i);
		// gmp_printf("(a,u,b,v)=(%Zd,%Zd,%Zd,%Zd)\n",a,u,b,v);
		mpz_sub(tempa, a, b); //sub a-b
		mpz_sub(tempu, u, v); //sub u-v
		// mpz_sub(tempa1, b, a); //sub b-a
		// mpz_sub(tempu1, v, u); //sub v-u
		s = mpz_tstbit(a, 0);
		result = mpz_cmp(a,b);// a>=bのときz=1 a<bのときz=0
		if (result >= 0) {
			z = 1;
		} else {
			z = 0;
		}

		b1 = (!s) || z;
		// printf("b1=%d\n",b1);
		//mpz_set(tempb, b); tempb->_mp_size *= b1;
		mpz_mul_ui(tempb, b, b1);
		//mpz_set(tempv, v); tempv->_mp_size *= b1;
		mpz_mul_ui(tempv, v, b1);
		// gmp_printf("tempb=%Zd\n",tempb);


		b1 = s && (!z);
		//printf("b1=%d\n",b1);
		//mpz_set(b,a);b->_mp_size *= b1;
		mpz_mul_ui(b, a, b1);
		//mpz_set(v,u);v->_mp_size *= b1;
		mpz_mul_ui(v, u, b1);
		mpz_xor(b, tempb, b);//xor
		mpz_xor(v, tempv, v);//xor
		mpz_mul_2exp(v, v, 1);//leftshift
		// tempa1->_mp_size *= b1;
		// result = (int)b1*(-1);
		// mpz_set_si(tempb, result);
		// gmp_printf("tempa=%Zd\n",tempa);
		// // gmp_printf("tempb=%Zd\n",tempb);
		// mpz_mul(tempa1, tempa, tempb);
		// gmp_printf("tempa1=%Zd\n",tempa1);
		// //tempu1->_mp_size *= b1;
		// mpz_mul(tempu1, tempu, tempb);

		result=~b1+1;
		mpz_set_si(tempb, result);
		// printf("b1=%d\n",b1);
		// printf("result=%d\n",result);
		// gmp_printf("tempa=%Zd\n",tempa);
		mpz_mul(tempa1, tempa, tempb);
		// gmp_printf("tempa1=%Zd\n",tempa1);
		mpz_mul(tempu1, tempu, tempb);

		b1 = !s;
		//printf("b1=%d\n",b1);
		//a->_mp_size *= b1;
		mpz_mul_ui(a, a, b1);
		mpz_xor(a, a, tempa1);
		//u->_mp_size *= b1;
		mpz_mul_ui(u, u, b1);
		mpz_xor(u, u, tempu1);
		
		b1 = s && z;
		//printf("b1=%d\n",b1);
		//tempa->_mp_size *= b1;
		mpz_mul_ui(tempa, tempa, b1);
		mpz_xor(a, a, tempa);
		// tempa xor a
		//tempu->_mp_size *= b1;
		mpz_mul_ui(tempu, tempu, b1);
		mpz_xor(u, u, tempu);
		// tempu xor u
		mpz_tdiv_q_2exp(a, a, 1);//rightshift
		
		// gmp_printf("(a,b)=(%Zd,%Zd)\n",a,b);
		//gmp_printf("(a,b,u,v)=(%Zd,%Zd,%Zd,%Zd)\n",a,b,u,v);

	}
	// gmp_printf("(u,v)=(%Zd,%Zd)\n",u,v);
	mpz_mul(v, v, pre_com);	mpz_mod(v, v, m);
	// gmp_printf("逆元=%Zd\n", v);

	mpz_clear(a); mpz_clear(b); mpz_clear(u); 
	mpz_clear(tempa); mpz_clear(tempa1); mpz_clear(tempb); 
	mpz_clear(tempu); mpz_clear(tempu1); mpz_clear(tempv); 
	return ;
}

void ex_test(mpz_t v, mpz_t y, mpz_t m, mpz_t pre_com,int LOOPS, mpz_t x1){
	mpz_t a,b,u,b2;
	mpz_init_set(a, y);
	mpz_init_set(b, m);
	mpz_init(b2);
	mpz_init_set_str(u, "1", 10);
	mpz_init_set_str(v, "0", 10);
	mpz_t tempa; mpz_init(tempa);
	mpz_t tempb; mpz_init(tempb);
	mpz_t tempu; mpz_init(tempu);
	mpz_t tempv; mpz_init(tempv);
	mpz_t tempa1; mpz_init(tempa1);
	mpz_t tempu1; mpz_init(tempu1);
	bool s, z, b1;
	int result;
	//gmp_printf("(a,u,b,v)=(%Zd,%Zd,%Zd,%Zd)\n",a,u,b,v);
	for (int i = 0; i < LOOPS; ++i){
		// printf("%d\n",i);
		// gmp_printf("(a,u,b,v)=(%Zd,%Zd,%Zd,%Zd)\n",a,u,b,v);
		mpz_sub(tempa, a, b); //sub a-b
		mpz_sub(tempu, u, v); //sub u-v
		// mpz_sub(tempa1, b, a); //sub b-a
		// mpz_sub(tempu1, v, u); //sub v-u
		s = mpz_tstbit(a, 0);
		result = mpz_cmp(a,b);// a>=bのときz=1 a<bのときz=0
		if (result >= 0) {
			z = 1;
		} else {
			z = 0;
		}

		b1 = (!s) || z;
		// printf("b1=%d\n",b1);
		if(b1 == 1){
			mpz_set(b2,x1);
		}else{
			mpz_set_si(b2,0);
		}
		// gmp_printf("b2=%Zd\n",b2);
		//printf("b1=%d\n",b1);
		//mpz_set(tempb, b); tempb->_mp_size *= b1;
		mpz_and(tempb, b, b2);
		//mpz_set(tempv, v); tempv->_mp_size *= b1;
		mpz_and(tempv, v, b2);
		// gmp_printf("tempb=%Zd\n",tempb);


		b1 = s && (!z);
		//printf("b1=%d\n",b1);
		if(b1 == 1){
			mpz_set(b2,x1);
		}else{
			mpz_set_si(b2,0);
		}
		//mpz_set(b,a);b->_mp_size *= b1;
		mpz_and(b, a, b2);
		//mpz_set(v,u);v->_mp_size *= b1;
		mpz_and(v, u, b2);
		mpz_xor(b, tempb, b);//xor
		mpz_xor(v, tempv, v);//xor
		mpz_mul_2exp(v, v, 1);//leftshift
		// tempa1->_mp_size *= b1;
		// result = (int)b1*(-1);
		// mpz_set_si(tempb, result);
		// gmp_printf("tempa=%Zd\n",tempa);
		// // gmp_printf("tempb=%Zd\n",tempb);
		// mpz_mul(tempa1, tempa, tempb);
		// gmp_printf("tempa1=%Zd\n",tempa1);
		// //tempu1->_mp_size *= b1;
		// mpz_mul(tempu1, tempu, tempb);

		result=~b1+1;
		mpz_set_si(tempb, result);
		// printf("b1=%d\n",b1);
		// printf("result=%d\n",result);
		// gmp_printf("tempa=%Zd\n",tempa);
		mpz_mul(tempa1, tempa, tempb);
		// gmp_printf("tempa1=%Zd\n",tempa1);
		mpz_mul(tempu1, tempu, tempb);

		b1 = !s;
		//printf("b1=%d\n",b1);
		if(b1 == 1){
			mpz_set(b2,x1);
		}else{
			mpz_set_si(b2,0);
		}
		//a->_mp_size *= b1;
		mpz_and(a, a, b2);
		mpz_xor(a, a, tempa1);
		//u->_mp_size *= b1;
		mpz_and(u, u, b2);
		mpz_xor(u, u, tempu1);
		
		b1 = s && z;
		//printf("b1=%d\n",b1);
		if(b1 == 1){
			mpz_set(b2,x1);
		}else{
			mpz_set_si(b2,0);
		}
		//tempa->_mp_size *= b1;
		mpz_and(tempa, tempa, b2);
		mpz_xor(a, a, tempa);
		// tempa xor a
		//tempu->_mp_size *= b1;
		mpz_and(tempu, tempu, b2);
		mpz_xor(u, u, tempu);
		// tempu xor u
		mpz_tdiv_q_2exp(a, a, 1);//rightshift
		
		// gmp_printf("(a,b)=(%Zd,%Zd)\n",a,b);
		//gmp_printf("(a,b,u,v)=(%Zd,%Zd,%Zd,%Zd)\n",a,b,u,v);

	}
	// gmp_printf("(u,v)=(%Zd,%Zd)\n",u,v);
	mpz_mul(v, v, pre_com);	mpz_mod(v, v, m);
	// gmp_printf("逆元=%Zd\n", v);

	mpz_clear(a); mpz_clear(b); mpz_clear(u); 
	mpz_clear(tempa); mpz_clear(tempa1); mpz_clear(tempb); 
	mpz_clear(tempu); mpz_clear(tempu1); mpz_clear(tempv); 
	return ;
}

void ex_test2(mpz_t v, mpz_t y, mpz_t m, mpz_t pre_com,int LOOPS, mpz_t x1){
	mpz_t a,b,u,b2;
	mpz_init_set(a, y);
	mpz_init_set(b, m);
	mpz_init_set_str(b2, "0", 10);
	mpz_init_set_str(u, "1", 10);
	mpz_init_set_str(v, "0", 10);
	mpz_t tempa; mpz_init(tempa);
	mpz_t tempb; mpz_init(tempb);
	mpz_t tempu; mpz_init(tempu);
	mpz_t tempv; mpz_init(tempv);
	mpz_t tempa1; mpz_init(tempa1);
	mpz_t tempu1; mpz_init(tempu1);
	bool s, z, b1;
	int result;
	//gmp_printf("(a,u,b,v)=(%Zd,%Zd,%Zd,%Zd)\n",a,u,b,v);
	for (int i = 0; i < LOOPS; ++i){
		// printf("%d\n",i);
		// gmp_printf("(a,u,b,v)=(%Zd,%Zd,%Zd,%Zd)\n",a,u,b,v);
		mpz_sub(tempa, a, b); //sub a-b
		mpz_sub(tempu, u, v); //sub u-v
		// mpz_sub(tempa1, b, a); //sub b-a
		// mpz_sub(tempu1, v, u); //sub v-u
		s = mpz_tstbit(a, 0);
		result = mpz_cmp(a,b);// a>=bのときz=1 a<bのときz=0
		if (result >= 0) {
			z = 1;
		} else {
			z = 0;
		}

		b1 = (!s) || z;
		// printf("b1=%d\n",b1);
		if(b1 == 1){
			mpz_ior(tempb, b, b2);
			mpz_ior(tempv, v, b2);
		}else{
			mpz_and(tempb,b ,b2);
			mpz_and(tempv,v ,b2);
		}
		// gmp_printf("b2=%Zd\n",b2);
		//printf("b1=%d\n",b1);
		//mpz_set(tempb, b); tempb->_mp_size *= b1;
		// mpz_and(tempb, b, b2);
		//mpz_set(tempv, v); tempv->_mp_size *= b1;
		// mpz_and(tempv, v, b2);
		// gmp_printf("tempb=%Zd\n",tempb);


		b1 = s && (!z);
		//printf("b1=%d\n",b1);
		if(b1 == 1){
			mpz_ior(b, a, b2);
			mpz_ior(v, u, b2);
		}else{
			mpz_and(b, a, b2);
			mpz_and(v, u, b2);
		}
		//mpz_set(b,a);b->_mp_size *= b1;
		// mpz_and(b, a, b2);
		//mpz_set(v,u);v->_mp_size *= b1;
		// mpz_and(v, u, b2);
		mpz_xor(b, tempb, b);//xor
		mpz_xor(v, tempv, v);//xor
		mpz_mul_2exp(v, v, 1);//leftshift
		// tempa1->_mp_size *= b1;
		// result = (int)b1*(-1);
		// mpz_set_si(tempb, result);
		// gmp_printf("tempa=%Zd\n",tempa);
		// // gmp_printf("tempb=%Zd\n",tempb);
		// mpz_mul(tempa1, tempa, tempb);
		// gmp_printf("tempa1=%Zd\n",tempa1);
		// //tempu1->_mp_size *= b1;
		// mpz_mul(tempu1, tempu, tempb);

		result=~b1+1;
		mpz_set_si(tempb, result);
		// printf("b1=%d\n",b1);
		// printf("result=%d\n",result);
		// gmp_printf("tempa=%Zd\n",tempa);
		mpz_mul(tempa1, tempa, tempb);
		// gmp_printf("tempa1=%Zd\n",tempa1);
		mpz_mul(tempu1, tempu, tempb);

		b1 = !s;
		//printf("b1=%d\n",b1);
		if(b1 == 1){
			mpz_ior(a, a, b2);
			mpz_ior(u, u, b2);
		}else{
			mpz_and(a, a, b2);
			mpz_and(u, u, b2);
		}
		//a->_mp_size *= b1;
		// mpz_and(a, a, b2);
		mpz_xor(a, a, tempa1);
		//u->_mp_size *= b1;
		// mpz_and(u, u, b2);
		mpz_xor(u, u, tempu1);
		
		b1 = s && z;
		//printf("b1=%d\n",b1);
		if(b1 == 1){
			mpz_ior(tempa, tempa, b2);
			mpz_ior(tempu, tempu, b2);
		}else{
			mpz_and(tempa, tempa, b2);
			mpz_and(tempu, tempu, b2);
		}
		//tempa->_mp_size *= b1;
		// mpz_and(tempa, tempa, b2);
		mpz_xor(a, a, tempa);
		// tempa xor a
		//tempu->_mp_size *= b1;
		// mpz_and(tempu, tempu, b2);
		mpz_xor(u, u, tempu);
		// tempu xor u
		mpz_tdiv_q_2exp(a, a, 1);//rightshift
		
		// gmp_printf("(a,b)=(%Zd,%Zd)\n",a,b);
		//gmp_printf("(a,b,u,v)=(%Zd,%Zd,%Zd,%Zd)\n",a,b,u,v);

	}
	// gmp_printf("(u,v)=(%Zd,%Zd)\n",u,v);
	mpz_mul(v, v, pre_com);	mpz_mod(v, v, m);
	// printf("ex_test2\n");
	gmp_printf("逆元=%Zd\n", v);

	mpz_clear(a); mpz_clear(b); mpz_clear(u); 
	mpz_clear(tempa); mpz_clear(tempa1); mpz_clear(tempb); 
	mpz_clear(tempu); mpz_clear(tempu1); mpz_clear(tempv); 
	return ;
}

void ex_test3(mpz_t *q, mpz_t a, mpz_t p, mpz_t pre_comp, int LOOPS){
	mpz_t u, v, r, temp1, temp2, temp3, temp4, temp5;
	mpz_init_set(u, a); mpz_init_set(v, p); mpz_set_ui(*q, 0); mpz_init_set_ui(r, 1);
	mpz_init(temp1); mpz_init(temp2); mpz_init(temp3); mpz_init(temp4); mpz_init(temp5);

	bool s, z;
	int result;
	// unsigned long int x;
	// x=1;
	
	mpz_t *A[3][2]; mpz_t *B[3][2];

	A[0][0] = A[2][1] = &u;
	A[0][1] = A[1][1] = &v;
	A[2][0] = &temp1;
	A[1][0] = &temp3;

	B[0][0] = B[2][1] = &r;
	B[0][1] = B[1][1] = q;
	B[2][0] = &temp2;
	B[1][0] = &temp4;
	
	
	uint8_t C[2][2];
	C[0][0] = C[0][1] = 0;
	C[1][0] = 2; C[1][1] = 1;

	for (int i = 0; i < LOOPS; ++i){
		// printf("%d回目\n",i);
		mpz_sub(temp1, v, u); //sub temp1 = v - u
		mpz_sub(temp2, *q, r); //sub temp2 = q - r
		mpz_com(temp3, temp1);
		mpz_add_ui(temp3, temp3, 1);
		// mpz_sub(temp3, u, v);
		mpz_com(temp4, temp2);
		mpz_add_ui(temp4, temp4, 1);
		// mpz_sub(temp4, r, *q);
		// gmp_printf("(u,v,q,r)=(%Zd,%Zd,%Zd,%Zd)\n",u,v,*q,r);
		// gmp_printf("(temp1,temp2)=(%Zd,%Zd)\n",temp1,temp2);
		// gmp_printf("(temp3,temp4)=(%Zd,%Zd)\n",temp3,temp4);
		s = mpz_tstbit(u, 0);
		result = mpz_cmp(u,v);// u>=vのときz=1 u<vのときz=0
		if (result >= 0) {
			z = 1;
		} else {
			z = 0;
		}
		
		C[0][0] = C[s][z]; 
		// printf("c[0][0]=%d\n",C[0][0]);

		mpz_tdiv_q_2exp(*A[C[0][0]][0], *A[C[0][0]][0], 1); //shift

		mpz_mul_2exp(*B[C[0][0]][1], *B[C[0][0]][1], 1); //shift

		// bs = mpz_cmp(*A[C[0][0]][0], *A[C[0][0]][1]) + 1; //32bit_add
		// bf = !bs;

	
		mpz_set(temp5, *A[C[0][0]][1]); // new v
		mpz_set(u, *A[C[0][0]][0]); // new u
		mpz_set(v, temp5);
		
		// mpz_set(v, temp5);
		// gmp_printf("(u,v,q,r)=(%Zd,%Zd,%Zd,%Zd)\n",u,v,*q,r);
		// mpz_set(temp5, *B[C[0][0]][0]); //new q
		// mpz_set(*q, temp5);
		// mpz_set(r, *B[C[0][0]][1]); //new r
		mpz_set(temp5, *B[C[0][0]][1]); //new q
		mpz_set(r, *B[C[0][0]][0]); // ner r
		mpz_set(*q, temp5);
		// gmp_printf("(u,v,q,r)=(%Zd,%Zd,%Zd,%Zd)\n",u,v,*q,r);
		
	}
	mpz_mul(*q, *q, pre_comp); mpz_mod(*q, *q, p); 
	// gmp_printf("逆元=%Zd\n", *q);
	mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3); 
	mpz_clear(u); mpz_clear(v); mpz_clear(r);
	mpz_clear(temp4); mpz_clear(temp5);
	return;
}


void ECTMI2(mpz_t q, mpz_t a, mpz_t p, mpz_t pre_comp, int LOOPS){
	mpz_t u, v, r;
	mpz_init_set(u, a); mpz_init_set(v, p); mpz_set_ui(q, 0); mpz_init_set_ui(r, 1);
	
	mpz_t temp1; mpz_init(temp1);
	mpz_t temp2; mpz_init(temp2);
	mpz_t temp3; mpz_init(temp3);
	bool u_lsb, v_lsb, t2;
	int t1;
	
	mpz_t *sort1[2]; sort1[0] = &temp1; sort1[1] = &temp2;
	mpz_t *sort2[2]; sort2[0] = &v; sort2[1] = &temp3;
	for (int i = 0; i < LOOPS; ++i){
		u_lsb = mpz_tstbit(u, 0);
		v_lsb = mpz_tstbit(v, 0);
		
		mpz_set(temp2, v); temp2->_mp_size *= u_lsb;
		//mpz_mul_ui(temp2, v, u_lsb);
		t1 = 2-(u_lsb<<1)-v_lsb;
		mpz_set(temp1, u); temp1->_mp_size *= t1;
		//mpz_mul_si(temp1, u, t1);
		mpz_add(temp2, temp2, temp1);// add
		mpz_tdiv_q_2exp(temp2, temp2, 1);// shift
		
		mpz_set(temp1, q); temp1->_mp_size *= u_lsb;
		//mpz_mul_ui(temp1, q, u_lsb);
		mpz_set(temp3, r); temp3->_mp_size *= t1;
		//mpz_mul_si(temp3, r, t1);
		mpz_add(temp3, temp1, temp3);//add
		
		t2 = v_lsb ^ u_lsb;
		v->_mp_size *= t2;
		//mpz_mul_ui(v, v, t2);
		t1 = ((v_lsb & u_lsb)<<1)-1;
		u->_mp_size *= t1;
		//mpz_mul_si(u, u, t1);
		mpz_add(temp1, v, u);//add
		
		q->_mp_size *= t2;
		//mpz_mul_ui(q, q, t2);
		r->_mp_size *= t1;
		//mpz_mul_si(r, r, t1);
		mpz_add(v, q, r);//add
		mpz_mul_2exp(v, v, 1);//shift

		u_lsb = mpz_cmp(temp2, temp1) + 1;
		t2 = !u_lsb;
		mpz_set(q, *sort2[u_lsb]);
		mpz_set(r, *sort2[t2]);
		mpz_set(v, *sort1[u_lsb]);
		mpz_set(u, *sort1[t2]);
		/*
		mpz_xor(u, temp1, temp2); //xor
		mpz_mul_ui(u, u, u_lsb); 
		mpz_xor(v, u, temp1); //xor
		mpz_xor(u, u, temp2); //xor
		*/
		/*
		mpz_mul_ui(v, temp2, b1);
		mpz_mul_ui(u, temp1, !b1);
		mpz_xor(v, v, u); //xor
		mpz_mul_ui(u, temp2, !b1);
		mpz_mul_ui(temp1, temp1, b1);
		mpz_xor(u, u, temp1); //xor
		*/
	}
	mpz_mul(q, q, pre_comp); mpz_mod(q, q, p); 
	mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3); 
	mpz_clear(u); mpz_clear(v); mpz_clear(r);
	return;
}


void KM1(mpz_t y, mpz_t a, mpz_t p, mpz_t pre_comp, int LOOPS){
	mpz_t u, v, x;
	mpz_init_set(u, a); mpz_init_set(v, p); mpz_set_ui(y, 0); mpz_init_set_ui(x, 1);
	
	mpz_t tempu; mpz_init(tempu);
	mpz_t tempv; mpz_init(tempv);
	mpz_t tempx; mpz_init(tempx);
	mpz_t tempy; mpz_init(tempy);

	bool s, z,  b1;
	int t1;
	int result;
	// gmp_printf("(u,v,x,y)=(%Zd,%Zd,%Zd,%Zd)\n",u,v,x,y);
	for (int i = 0; i < LOOPS; ++i){
		// printf("%d回目\n",i);
		s = mpz_tstbit(u, 0);
		// v_lsb = mpz_tstbit(v, 0);
		result = mpz_cmp(u,v);// u>=vのときz=1 u<vのときz=0
		if (result >= 0) {
			z = 1;
		} else {
			z = 0;
		}
		// printf("s,z=%d,%d\n",s,z);
		b1 = (!s) || z;
		// printf("b1=%d\n",b1);
		// u,xを求めるフェーズ
		t1 = 2-(b1<<1)-s;
		// printf("b1=%d\n",b1<<1);
		// printf("s,z,sz=%d,%d,%d\n",s,z,s&z);
		// printf("t1=%d\n",t1);
		mpz_mul_si(tempv, v, t1);
		mpz_mul_si(tempy, y, t1);
		t1 = (b1<<1) - 1;
		// printf("t1=%d\n",t1);
		mpz_mul_si(tempu, u, t1);
		mpz_mul_si(tempx, x, t1);
		// gmp_printf("tempu=%Zd,",tempu);
		mpz_add(tempu, tempv, tempu);// add
		mpz_add(tempx, tempx, tempy);// add temp_x =new x
	
		mpz_tdiv_q_2exp(tempu, tempu, 1);// shift temp_u =new u


		// v,yを求めるフェーズ
		mpz_mul_ui(v, v, b1);
		// gmp_printf("v=%Zd,",v);
		mpz_mul_ui(y, y, b1);

		b1 = !b1;
		mpz_mul_ui(u, u, b1);
		// gmp_printf("u=%Zd,",u);
		mpz_mul_ui(x, x, b1);
		mpz_add(v, u, v);//add
		mpz_add(y, x, y);//add
		mpz_mul_2exp(y, y, 1);//shift

		mpz_set(u, tempu);
		mpz_set(x, tempx);
		// gmp_printf("(u,v,x,y)=(%Zd,%Zd,%Zd,%Zd)\n",u,v,x,y);

	}
	mpz_mul(y, y, pre_comp); mpz_mod(y, y, p); 
	gmp_printf("逆元=%Zd\n", y);
	mpz_clear(tempu); mpz_clear(tempv); mpz_clear(tempx); mpz_clear(tempy);
	mpz_clear(u); mpz_clear(v); mpz_clear(x);
	return;
}


void ECTMI3(mpz_t *q, mpz_t a, mpz_t p, mpz_t pre_comp, int LOOPS){
	mpz_t u, v, r, temp1, temp2;
	mpz_init_set(u, a); mpz_init_set(v, p); mpz_set_ui(*q, 0); mpz_init_set_ui(r, 1);
	mpz_init(temp1); mpz_init(temp2);

	bool u_lsb, v_lsb, z;
	int result;
	mpz_t temp3; mpz_init(temp3);
	
	mpz_t *A[3][2]; mpz_t *B[3][2];

	A[0][0] = A[1][1] = A[2][1] = &temp1;
	A[1][0] = &v;
	A[2][0] = A[0][1] = &u;
	
	B[0][0] = B[1][1] = B[2][1] = &temp2;
	B[1][0] = q;
	B[2][0] = B[0][1] = &r;
	
	uint8_t C[2][2];
	C[0][0] = C[1][1] = 0;
	C[0][1] = 2; C[1][0] = 1;

	for (int i = 0; i < LOOPS; ++i){
		mpz_sub(temp1, v, u); //sub temp1 = v - u
		mpz_sub(temp2, *q, r); //sub temp2 = q - r

		u_lsb = mpz_tstbit(u, 0);
		v_lsb = mpz_tstbit(v, 0);
		
		C[0][0] = C[u_lsb][v_lsb]; 

		mpz_tdiv_q_2exp(*A[C[0][0]][0], *A[C[0][0]][0], 1); //shift

		mpz_mul_2exp(*B[C[0][0]][1], *B[C[0][0]][1], 1); //shift

		result = mpz_cmp(*A[C[0][0]][0], *A[C[0][0]][1]); //mpz_cmp(u,v)u>=vのときz=0 u<vのときz=1
		if (result >= 0) {
			z = 0;
		} else {
			z = 1;
		}

		mpz_set(temp3, *A[C[0][0]][z]); // new v
		mpz_set(u, *A[C[0][0]][!z]); 
		mpz_set(v, temp3);

		mpz_set(temp3, *B[C[0][0]][z]); //new q
		mpz_set(r, *B[C[0][0]][!z]);
		mpz_set(*q, temp3);
	}
	mpz_mul(*q, *q, pre_comp); mpz_mod(*q, *q, p); 
	// gmp_printf("逆元=%Zd\n", *q);
	mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3); 
	mpz_clear(u); mpz_clear(v); mpz_clear(r);
	return;
}


// void BOS(mpz_t *q, mpz_t a, mpz_t p, int LOOPS){
// 	mpz_t z,u,v,r,k,d;
// 	mpz_t zero,one;
// 	mpz_t m1,m2,m3,m4,m5,m6,m7,S;
// 	mpz_init(m1); mpz_init(m2); mpz_init(m3); mpz_init(m4); mpz_init(m5); mpz_init(m6); mpz_init(m7);
// 	mpz_init(S);
// 	mpz_init_set(u, a);
// 	mpz_init_set(v, p);
// 	mpz_init_set_str(r, "1", 10);
// 	mpz_init_set_str(q, "0", 10);
// 	mpz_init_set_str(k, "0",10);
// 	mpz_init_set_str(zero, "0",10);
// 	mpz_init_set_str(one, "1",10);

// 	mpz_t sub_vu; mpz_init(sub_vu);
// 	int vu;


// 	bool u_lsb,v_lsb;
	
// 	//gmp_printf("(a,u,b,v)=(%Zd,%Zd,%Zd,%Zd)\n",a,u,b,v);
// 	for (int i = 0; i < LOOPS; ++i){
// 		u_lsb = mpz_tstbit(u, 0); //uの最下位ビット
// 		v_lsb = mpz_tstbit(v, 0); //vの最下位ビット
// 		mpz_sub(sub_vu, v, u); //sub v-u
// 		vu = mpz_cmp(v,u);// v>=uのときz=1 a<bのときz=0
// 		if (vu > 0) {
// 			mpz_set(z,1);
// 		} else {
// 			mpz_set(z,0);
// 		}
// 		mpz_sub_uni(d, 0, sub_vu);
// 		mpz_tdiv_q_2exp(v, v, 1); //v=v>>1
// 		mpz_mul_2exp(r, r, 1); //r=r<<1
// 		mpz_and(u_lsb, u_lsb, one); // u_lsb and 1
// 		mpz_sub(u_lsb, 0, u_lsb); // u_lsb = 0 - u_lsb
// 		mpz_ior(m1);


// 	}
// 	// // gmp_printf("(u,v)=(%Zd,%Zd)\n",u,v);
// 	// mpz_mul(v, v, pre_com);	mpz_mod(v, v, m);
// 	// // gmp_printf("逆元=%Zd\n", v);

// 	// mpz_clear(a); mpz_clear(b); mpz_clear(u); 
// 	// mpz_clear(tempa); mpz_clear(tempa1); mpz_clear(tempb); 
// 	// mpz_clear(tempu); mpz_clear(tempu1); mpz_clear(tempv); 
// 	return ;
// }

static __inline__ uint64_t GetCC(){
	unsigned int a, d;
	asm volatile("mfence;rdtsc" : "=a" (a), "=d" (d)::"memory");
	return ((uint64_t)a) | (((uint64_t)d)<<32);
}

