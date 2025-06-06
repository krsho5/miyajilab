#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <gmp.h>

/*FCTMI: The number of loops defined by (49d+57)//17 when d = 224, 256, or 384.*/
#define P_224_LOOPS_FCTMI 634//649 
#define P_256_LOOPS_FCTMI 724//741
#define P_384_LOOPS_FCTMI 1086//1110

/*Precomputed numbers of FCTMI.*/
#define P_224_PRE_FCIMI "2008672555230201739572940614232879648470511825043254648242560"//"26644009792206268056205042062915547495932996820470046226332801564673"
#define P_256_PRE_FCIMI "115789218083876703796405708104922954266844366797366597945313458537426472636415"//"90462569673685862684289949614566734656231367827980882643939160353509669339137"
#define P_384_PRE_FCIMI "39401102631331629224042747884736582634501236456156400192182914717222615670170491983949444240325277863502310224573143"//"39124746653188674882124692525769249793887199085144063008056343542510001139863950072793525947872276405658285025530879"

/*hdFCTMI: The number of loops when d = 224, 256, or 384.*/
#define P_224_LOOPS_hdFCTMI 517 
#define P_256_LOOPS_hdFCTMI 590
#define P_384_LOOPS_hdFCTMI 885

/*ECTMI: The number of loops defined by 2*bitlen(p) with p224, p256, p384.*/
#define P_224_LOOPS_ECTMI 448 
#define P_256_LOOPS_ECTMI 512
#define P_384_LOOPS_ECTMI 768

/*The number of pre_com 2^(-l) mod p*/
#define P_256_pre_com  "161759680002903838760694582340907787320905632487091026015073739997175"

/*Definition of the field.*/
#define P_224 "26959946667150639794667015087019630673557916260026308143510066298881"
#define P_256 "115792089210356248762697446949407573530086143415290314195533631308867097853951"
#define P_384 "39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319"

/*BY19	GCD(a, p) with a or p is odd*/
void FCTMI(mpz_t v, mpz_t a, mpz_t p, int LOOPS){
	mpz_t u;
	mpz_init_set(u, a); mpz_set(v, p);

	int delta = 1; 
	int s, z; 
	int temp1; mpz_t temp2, temp3; mpz_init(temp2); mpz_init(temp3);

	for (int i = 0; i < LOOPS; ++i){
		z = mpz_tstbit(u, 0); //z = LSB(u)
		s = (delta >> 31) + 1; //s = signbit(-delta)

		/*delta = 1 + (1 - 2sz) delta*/
		s = s & z; // s = sz
		temp1 = (s << 1); // temp1 = 2sz
		temp1 = 1 - temp1; // temp1 = 1-2sz
		delta *= temp1; // delta = (1-2sz) delta
		delta += 1; // delta = delta + 1

		/*part 1 u = (u + (1-2sz)z v)/2*/
		z *= temp1; // z = (1-2sz)z
		mpz_set(temp2, v); temp2->_mp_size *= z;
		//mpz_mul_si(temp2, v, z); // temp2 = (1-2sz)zv
		
		/*v = v xor sz(v xor u)*/
		mpz_xor(temp3, v, u); // temp3 = v xor u
		temp3->_mp_size*=s;
		//mpz_mul_ui(temp3, temp3, s); // temp3 = sz(v xor u)
		mpz_xor(v, v, temp3); // v = v xor sz(v xor u)

		/*part 2 u = (u + (1-2sz)z v)/2*/
		mpz_add(u, u, temp2); // u = u + (1-2sz)zv
		mpz_tdiv_q_2exp(u, u, 1); // u = (u + (1-2sz)zv)/2
	}
	mpz_abs(v, v); 
	mpz_clear(u); mpz_clear(temp2); mpz_clear(temp3);
	return;
}

/*advanced BY19 GCD(a, p) with a or p is odd*/
void hdFCTMI(mpz_t v, mpz_t a, mpz_t p, int LOOPS){
	mpz_t u;
	mpz_init_set(u, a); mpz_set(v, p);

	int delta = 1; 
	int s, z; 
	int temp1; mpz_t temp2, temp3; mpz_init(temp2); mpz_init(temp3);

	for (int i = 0; i < LOOPS; ++i){
		z = mpz_tstbit(u, 0); //z = LSB(u)
		s = (delta >> 31) + 1; //s = signbit(-delta)

		/*delta = 2 + (1 - 2sz) delta*/
		s = s & z; // s = sz
		temp1 = (s << 1); // temp1 = 2sz
		temp1 = 1 - temp1; // temp1 = 1-2sz
		delta *= temp1; // delta = (1-2sz) delta
		delta += 2; // delta = delta + 2

		/*part 1 u = (u + (1-2sz)z v)/2*/
		z *= temp1; // z = (1-2sz)z
		mpz_set(temp2, v); temp2->_mp_size *= z;
		//mpz_mul_si(temp2, v, z); // temp2 = (1-2sz)zv
		
		/*v = v xor sz(v xor u)*/
		mpz_xor(temp3, v, u); // temp3 = v xor u
		temp3->_mp_size*=s;
		//mpz_mul_si(temp3, temp3, s); // temp3 = sz(v xor u)
		mpz_xor(v, v, temp3); // v = v xor sz(v xor u)

		/*part 2 u = (u + (1-2sz)z v)/2*/
		mpz_add(u, u, temp2); // u = u + (1-2sz)zv
		mpz_tdiv_q_2exp(u, u, 1); // u = (u + (1-2sz)zv)/2
	}
	mpz_abs(v, v); 
	mpz_clear(u); mpz_clear(temp2); mpz_clear(temp3);
	return;
}

/*Our GCD(a, p) with a or p is odd and a <= p*/
/*
void ECTMI_IF(mpz_t v, mpz_t a, mpz_t p, int LOOPS){
	mpz_t u;
	mpz_init_set(u, a); mpz_set(v, p);
	
	int u_lsb, v_lsb;
	mpz_t temp1, temp2; mpz_init(temp1); mpz_init(temp2);
	
	for (int i = 0; i < LOOPS; ++i){
		u_lsb = mpz_tstbit(u, 0);
		v_lsb = mpz_tstbit(v, 0);
		mpz_sub(temp1, v, u); //sub temp1 = v - u 
		if(u_lsb == 1 && v_lsb == 1){
			mpz_tdiv_q_2exp(temp1, temp1, 1); //shift
			if(mpz_cmp(u, temp1) > 0){
				mpz_set(temp2, u); mpz_set(v, temp2); mpz_set(temp2, temp1); mpz_set(u, temp2);// 4assign
			}
			else{
				mpz_set(temp2, u); mpz_set(u, temp1); mpz_set(v, u); mpz_set(u, temp2);// 4assign
			}
		}
		else if(u_lsb == 1 && v_lsb == 0){
			mpz_tdiv_q_2exp(u, v, 1); //shift
			if(mpz_cmp(u, temp1) > 0){
				mpz_set(temp2, u); mpz_set(v, temp2); mpz_set(temp2, temp1); mpz_set(u, temp2); // 4assign
			}
			else{
				mpz_set(temp2, u); mpz_set(u, temp1); mpz_set(v, u); mpz_set(u, temp2);// 4assign
			}
		}
		else{
			mpz_tdiv_q_2exp(u, u, 1); //shift
			if(mpz_cmp(u, temp1) > 0){
				mpz_set(temp2, u); mpz_set(v, temp2); mpz_set(temp2, temp1); mpz_set(u, temp2); // 4assign
			}
			else{
				mpz_set(temp2, u); mpz_set(u, temp1); mpz_set(v, u); mpz_set(u, temp2);// 4assign
			}
		}
	}
	mpz_clear(u); mpz_clear(temp1); mpz_clear(temp2);
	return;
}
*/

void ECTMI(mpz_t v, mpz_t a, mpz_t p, int LOOPS){
	mpz_t u;
	mpz_init_set(u, a); mpz_set(v, p);
	
	mpz_t temp1; mpz_init(temp1);
	mpz_t temp2; mpz_init(temp2);
	bool u_lsb, v_lsb, b1;
	
	mpz_t *sort[2]; sort[0] = &temp1; sort[1] = &temp2;
	for (int i = 0; i < LOOPS; ++i){
		mpz_sub(temp1, v, u); //sub

		u_lsb = mpz_tstbit(u, 0);
		v_lsb = mpz_tstbit(v, 0);
		
		b1 = (!u_lsb) && v_lsb;
		mpz_set(temp2, u); temp2->_mp_size *= b1;
		//mpz_mul_ui(temp2, u, b1);

		b1 = u_lsb && (!v_lsb);
		v->_mp_size *= b1;
		//mpz_mul_ui(v, v, b1);
		mpz_xor(temp2, temp2, v);//xor

		b1 = u_lsb && v_lsb;
		mpz_set(v, temp1); v->_mp_size *= b1;
		//mpz_mul_ui(v, temp1, b1);
		mpz_xor(temp2, temp2, v);//xor
		mpz_tdiv_q_2exp(temp2, temp2, 1);//shift
		u->_mp_size *= b1;
		//mpz_mul_ui(u, u, b1);

		b1 = (!u_lsb) || (!v_lsb);
		temp1->_mp_size *= b1;
		//mpz_mul_ui(temp1, temp1, b1);
		mpz_xor(temp1, temp1, u);//xor

		b1 = mpz_cmp(temp2, temp1) + 1;
		mpz_set(v, *sort[b1]);
		mpz_set(u, *sort[!b1]);
		/*
		mpz_xor(u, temp1, temp2); //xor
		mpz_mul_ui(u, u, b1); 
		mpz_xor(v, u, temp1); //xor
		mpz_xor(u, u, temp2); //xor
		*/
	}
	mpz_clear(u); mpz_clear(temp1); mpz_clear(temp2); 
	return;
}

void ECTMI2(mpz_t v, mpz_t a, mpz_t p, int LOOPS){
	mpz_t u;
	mpz_init_set(u, a); mpz_set(v, p);
	
	mpz_t temp1; mpz_init(temp1);
	mpz_t temp2; mpz_init(temp2);
	bool u_lsb, v_lsb, t2;
	int t1;
	
	mpz_t *sort[2]; sort[0] = &temp1; sort[1] = &temp2;
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
		
		t2 = v_lsb ^ u_lsb;
		v->_mp_size *= t2;
		//mpz_mul_ui(v, v, t2);
		t1 = ((v_lsb & u_lsb)<<1)-1;
		u->_mp_size *= t1;
		//mpz_mul_si(u, u, t1);
		mpz_add(temp1, v, u);//add

		u_lsb = mpz_cmp(temp2, temp1) + 1;
		mpz_set(v, *sort[u_lsb]);
		mpz_set(u, *sort[!u_lsb]);
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
	mpz_clear(u); mpz_clear(temp1); mpz_clear(temp2); 
	return;
}

void ECTMI3(mpz_t *v, mpz_t a, mpz_t p, int LOOPS){
	mpz_t u, temp1;
	mpz_init_set(u, a); mpz_set(*v, p);
	mpz_init(temp1);

	bool u_lsb, v_lsb, bf, bs;
	mpz_t temp3; mpz_init(temp3);
	
	mpz_t *A[3][2];

	A[0][0] = A[1][1] = A[2][1] = &temp1;
	A[1][0] = v;
	A[2][0] = A[0][1] = &u;
	
	uint8_t C[2][2];
	C[0][0] = C[1][1] = 0;
	C[0][1] = 2; C[1][0] = 1;

	for (int i = 0; i < LOOPS; ++i){
		mpz_sub(temp1, *v, u); //sub temp1 = v - u

		u_lsb = mpz_tstbit(u, 0);
		v_lsb = mpz_tstbit(*v, 0);
		
		C[0][0] = C[u_lsb][v_lsb]; 

		mpz_tdiv_q_2exp(*A[C[0][0]][0], *A[C[0][0]][0], 1); //shift

		bs = mpz_cmp(*A[C[0][0]][0], *A[C[0][0]][1]) + 1; //32bit_add
		bf = !bs;

		mpz_set(temp3, *A[C[0][0]][bf]); // new v
		mpz_set(u, *A[C[0][0]][bs]); 
		mpz_set(*v, temp3);
	}
	mpz_clear(u); mpz_clear(temp1); mpz_clear(temp3); 
	return;
}

static __inline__ uint64_t GetCC(){
	unsigned int a, d;
	asm volatile("mfence;rdtsc" : "=a" (a), "=d" (d)::"memory");
	return ((uint64_t)a) | (((uint64_t)d)<<32);
}

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
	for (int i = 0; i < LOOPS; ++i){
		mpz_sub(temp1, v, u); //sub v-u
		mpz_sub(temp3, q, r); //sub q-r
		u_lsb = mpz_tstbit(u, 0);
		v_lsb = mpz_tstbit(v, 0);
		
		b1 = (!u_lsb) && v_lsb;
		mpz_set(temp2, u); temp2->_mp_size *= b1;
		//mpz_mul_ui(temp2, u, b1);
		mpz_set(temp4, r); temp4->_mp_size *= b1;
		//mpz_mul_ui(temp4, r, b1);

		b1 = u_lsb && (!v_lsb);
		v->_mp_size *= b1;
		//mpz_mul_ui(v, v, b1);
		mpz_xor(temp2, temp2, v);//xor
		q->_mp_size *= b1;
		//mpz_mul_ui(q, q, b1);
		mpz_xor(temp4, temp4, q);//xor

		b1 = u_lsb && v_lsb;
		mpz_set(v, temp1); v->_mp_size *= b1;
		//mpz_mul_ui(v, temp1, b1);
		mpz_xor(temp2, temp2, v);//xor
		mpz_tdiv_q_2exp(temp2, temp2, 1);//rightshift
		u->_mp_size *= b1;
		//mpz_mul_ui(u, u, b1);
		mpz_set(q, temp3); q->_mp_size *= b1;
		//mpz_mul_ui(q, temp3, b1);
		mpz_xor(temp4, temp4, q);//xor
		r->_mp_size *= b1;
		//mpz_mul_ui(r, r, b1);

		b1 = (!u_lsb) || (!v_lsb);
		temp1->_mp_size *= b1;
		//mpz_mul_ui(temp1, temp1, b1);
		mpz_xor(temp1, temp1, u);//xor
		temp3->_mp_size *= b1;
		//mpz_mul_ui(temp3, temp3, b1);
		mpz_xor(temp3, temp3, r);//xor
		mpz_mul_2exp(temp3, temp3, 1);//leftshift
		
		b1 = mpz_cmp(temp2, temp1) + 1;
		mpz_set(v, *sort1[b1]); mpz_set(u, *sort1[!b1]);
		mpz_set(q, *sort2[b1]); mpz_set(r, *sort2[!b1]);
	}
	mpz_mul(q, q, pre_com);	mpz_mod(q, q, p);
	// gmp_printf("逆元=%Zd\n", q);
	mpz_clear(u); mpz_clear(temp1); mpz_clear(temp2); 
	mpz_clear(q); mpz_clear(r); mpz_clear(temp3); mpz_clear(temp4); 
	return;
}

void BY(mpz_t q, mpz_t a, mpz_t p, mpz_t pre_comp, int LOOPS){
	mpz_t u, v;
	mpz_init_set(u, a); mpz_init_set(v, p);

	mpz_t r; 
	mpz_set_ui(q, 0); mpz_init_set_ui(r, 1);

	int delta = 1; 
	int s, z; 
	int temp1; mpz_t temp2, temp3; mpz_init(temp2); mpz_init(temp3);

	for (int i = 0; i < LOOPS; ++i){
		z = mpz_tstbit(u, 0); //z = LSB(u)
		s = (delta >> 31) + 1; //s = signbit(-delta)

		/*delta = 1 + (1 - 2sz) delta*/
		s = s & z; // s = sz
		temp1 = (s << 1); // temp1 = 2sz
		temp1 = 1 - temp1; // temp1 = 1-2sz
		delta *= temp1; // delta = (1-2sz) delta
		delta += 1; // delta = delta + 1

		/*part 1 u = (u + (1-2sz)z v)/2*/
		z *= temp1; // z = (1-2sz)z
		// mpz_set(temp2, v); temp2->_mp_size *= z;
		mpz_mul_si(temp2, v, z); // temp2 = (1-2sz)zv
		
		/*v = v xor sz(v xor u)*/
		mpz_xor(temp3, v, u); // temp3 = v xor u
		// temp3->_mp_size *= s;
		mpz_mul_ui(temp3, temp3, s); // temp3 = sz(v xor u)
		mpz_xor(v, v, temp3); // v = v xor sz(v xor u)

		/*part 2 u = (u + (1-2sz)z v)/2*/
		mpz_add(u, u, temp2); // u = u + (1-2sz)zv
		mpz_tdiv_q_2exp(u, u, 1); // u = (u + (1-2sz)zv)/2

		/*part 1 r = (1-2sz) z q + r*/
		// mpz_set(temp2, q); temp2->_mp_size *= z;
		mpz_mul_si(temp2, q, z); // temp2 = (1-2sz)z q

		/*q = 2(q xor sz(q xor r)) q = 2(sz r+(sz XOR 1)q)*/
		mpz_xor(temp3, q, r); // temp3 = q xor r
		// temp3->_mp_size *= s;
		mpz_mul_ui(temp3, temp3, s); // temp3 = sz(q xor r)
		mpz_xor(q, q, temp3); // q = q xor sz(q xor r)
		mpz_mul_2exp(q, q, 1); // q = 2(q xor sz(q xor r)) 
		/*
		temp1 = s ^ 1; // temp1 = sz XOR 1
		mpz_mul_si(q, q, temp1); // q = q (sz XOR 1)
		mpz_mul_si(temp3, r, s); // temp3 = r sz
		mpz_add(q, q, temp3); // q = q (sz XOR 1) + r sz
		mpz_mul_2exp(q, q, 1); // q = 2(q (sz XOR 1) + r sz)
		*/

		/*part 2 r = (1-2sz) z q + r*/
		mpz_add(r, r, temp2);		
	}
	// q->_mp_size *= mpz_sgn(v);
	mpz_mul_si(q, q, mpz_sgn(v)); 
	mpz_mul(q, q, pre_comp); mpz_mod(q, q, p); 
	gmp_printf("逆元=%Zd\n", q);
	mpz_clear(u); mpz_clear(v); mpz_clear(r); mpz_clear(temp2); mpz_clear(temp3);
	return;
}
