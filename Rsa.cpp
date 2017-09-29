#include<iostream>
#include<gmpxx.h>
#include<ctime>
#include<cstdlib>
#include<unistd.h>

using namespace std;

/*
 * Secure Outsource RSA Decryption
 * Five Steps
 * 1. Preprocess
 * 2. Problem Transformation
 * 3. Compute
 * 4. Verify
 * 5. Recovery
 *
 */



// 全局变量
int randMpzFlag = 0;         // 标记randMpz()首次运行
int randFlag = 0;            // 标记random()首次运行
gmp_randstate_t rstate;      // 随机数生成参数rstate,仅初始化一次



/*
 * # Function
 * random()    : produce random number(int) in range of min~max
 * # Paramater
 * 
 */
int random(int min,int max){

    int ret;

    if(randFlag == 0){
        // initialize random seed:
        srand(time(NULL));
    }

    // generate secret number between min and max:
    ret = rand()%(max-min+1)+min;

    // mark the first run of random()
    randFlag = 1;

    return ret;
}



/*
 * # Function
 * randMpz()  : 生成并返回指定bit的随机数
 * # Parameter
 * retMpz     : 生成的随机大数
 * randBit    : 二进制位数
 * 注: gmp lib 中mpz_t不能作为函数返回值,gmp会自动回传参数retMPz
 */
void randMpz(mpz_t retMpz,unsigned long randBit){

    unsigned long seed;    // 初始化seed,调用time(NULL)根据当前时间生成seed
    if(randMpzFlag==0){    // 若首次运行则初始化 seed 及 rstate,若每次初始化,那么随机数大小相近
        seed = time(NULL);
        gmp_randinit_default(rstate);
        gmp_randseed_ui(rstate,seed);
    }

    mpz_t start;
    mpz_t end;
    mpz_t tmp;

    mpz_init(start);
    mpz_init(end);
    mpz_init(tmp);

    // 求随机数区间start~end
    mpz_ui_pow_ui(start,2,randBit-1);
    mpz_ui_pow_ui(end,2,randBit);
    mpz_sub_ui(end,end,1);

    // 求随机数
    mpz_urandomb(retMpz,rstate,randBit);
    mpz_sub(tmp,end,start);
    mpz_add_ui(tmp,tmp,1);
    mpz_mod(retMpz,retMpz,tmp);
    mpz_add(retMpz,retMpz,start);

    // 标记randFlag,下次不再初始化 seed 和 rstate
    randMpzFlag=1;

    // 清空内存
    mpz_clear(start);
    mpz_clear(end);
    mpz_clear(tmp);
}



/*
 * # Function
 * primeMpz   : 生成大素数
 * # Parameter
 * rand       : 传入的随机大数，以此为基准寻找附近的大素数,随后返回值也是rand
 * Miller-Rabin正确率不依赖被检测数p，而仅依赖检测次数k
 * k次检测后得到错误结果的概率为(1/4)^k,这里取k=50
 * Reasonable values of reps are between 15 and 50
 */
void primeMpz(mpz_t rand){

    // 如果rand为偶数,则加1变为奇数
    if(mpz_even_p(rand)){
        mpz_add_ui(rand,rand,1);
    }
    while(mpz_probab_prime_p(rand,50)==0){
        mpz_add_ui(rand,rand,2);
    }
}



/*
 * 计算所用时间
 */
double timeCost(timespec start, timespec end){

    timespec temp;
    if ((end.tv_nsec-start.tv_nsec)<0) {
        temp.tv_sec = end.tv_sec-start.tv_sec-1;
        temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
        temp.tv_sec = end.tv_sec-start.tv_sec;
        temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    double ret;
    ret = (double)temp.tv_sec + (double)temp.tv_nsec/1000000000;
    return ret;
}



/*
 * Entrance of main function
 */
int main(){

    int round = 0;                 // The number of cycles.
    int count = 0;                 // Record the round.

    int cmp1;                      // Record Verify result.
    int cmp2;

    int state;                     // Record the result state of mpz_invert()

    double avgTimeA1 = 0.0;        // Record the time cost of each phase.
    double avgTimeA2 = 0.0;
    double avgTimeA3 = 0.0;
    double avgTimeA4 = 0.0;
    double avgTimeA5 = 0.0;
    double avgTimeA6 = 0.0;

    double avgTimeB1 = 0.0;        // Record the time cost of each phase.
    double avgTimeB2 = 0.0;
    double avgTimeB3 = 0.0;
    double avgTimeB4 = 0.0;
    double avgTimeB5 = 0.0;
    double avgTimeB6 = 0.0;



    /*
     * Define the bit length of variables.
     */
    unsigned long pLen;           // The length of p and q are 1/2 of C.
    unsigned long qLen;
    
    unsigned long CLen;           // User input the length of C and d.
    unsigned long dLen;
    
    unsigned long rLen = 80;      // The length of r1,r2,r3 and r4 are fixed.
    unsigned long kLen = 80;

    unsigned long t1;             // t1,t2 and k are in the range of 2~11.
    unsigned long t2;
    unsigned long k;

    unsigned long r,t,w;          // Used in another plan.


    timespec beginT,endT;         // Record the time point.



    /*
     * Define the numbers of mpz_t.
     */
    mpz_t L;
    mpz_t k1,k2;                  // k1,k2 : 80bits
    mpz_t d,dd1,dd2;              // dd1 = d1' , dd2 = d2'
    mpz_t M,C,CC;                 // CC  = C'
    mpz_t p,q,n;                  // n   = pq.
    mpz_t p1,q1;                  // p1  = p'  ,  q1 = q'
    mpz_t dp,dq;
    mpz_t Mp,Mq;
    mpz_t Y1,Y2;
    mpz_t p_1,q_1,n_1;            // Record p-1,q-1,n-1.
    mpz_t r1,r2,r3,r4;
    mpz_t dp1,dp2,dq1,dq2;
    mpz_t R1,RR1,R2,RR2;          // RR1 = R1' , RR2 = R2'
    mpz_t tmp1,tmp2,tmp3;
    mpz_t tmp4,tmp5,tmp6;



    /*
     * Initial the number of mpz_t.
     */
    mpz_init(L);

    mpz_init(k1);
    mpz_init(k2);

    mpz_init(d);
    mpz_init(dd1);
    mpz_init(dd2);

    mpz_init(M);
    mpz_init(C);
    mpz_init(CC);

    mpz_init(p);
    mpz_init(q);
    mpz_init(n);

    mpz_init(p1);
    mpz_init(q1);

    mpz_init(dp);
    mpz_init(dq);

    mpz_init(Mp);
    mpz_init(Mq);

    mpz_init(Y1);
    mpz_init(Y2);

    mpz_init(p_1);
    mpz_init(q_1);
    mpz_init(n_1);

    //mpz_init(r);
    mpz_init(r1);
    mpz_init(r2);
    mpz_init(r3);
    mpz_init(r4);

    mpz_init(dp1);
    mpz_init(dp2);
    mpz_init(dq1);
    mpz_init(dq2);

    mpz_init(R1);
    mpz_init(RR1);
    mpz_init(R2);
    mpz_init(RR2);

    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init(tmp3);
    mpz_init(tmp4);
    mpz_init(tmp5);
    mpz_init(tmp6);



    /*
     * User input bit length of mpz_t.
     */
    cout<<endl;
    cout<<"Default: n = pq                                  "<<endl;
    cout<<"         L = pn                                  "<<endl;
    cout<<"         p and q                     : 1/2  CLen."<<endl;
    cout<<"         r1,r2,r3 and r4             : 80   bits."<<endl;
    cout<<"         k1 and k2                   : 80   bits."<<endl;
    cout<<"         k,t1,t2,r,t,w range in 2~11.      "<<endl<<endl;
    cout<<"Input the bit length of C (256~2048) :";
    cin>>CLen;
    cout<<"Input the bit length of d (256~2048) :";
    cin>>dLen;
    cout<<"Input the round of cycles (Integer ) :";
    cin>>round;
    cout<<endl;

    pLen = CLen/2;
    qLen = CLen/2;



    // Show the length of the progress bar.
    for( count=0; count<round; count++ ){
        printf(".");
    }
    printf("\n");



    for( count=0; count<round; count++ ){
        // Show the progress bar.
        cout<<"."<<flush;

        /*
         * Invoke randMpz() to produce random number
         */
        randMpz(d,dLen);
        randMpz(p,pLen);
        randMpz(q,qLen);

        randMpz(C,CLen);
        randMpz(d,dLen);

        randMpz(k1,kLen);
        randMpz(k2,kLen);

        //randMpz(r,rLen);
        randMpz(r1,rLen);
        randMpz(r2,rLen);
        randMpz(r3,rLen);
        randMpz(r4,rLen);



        /*
         * Produce a Prime based on the random p,q
         */
        primeMpz(p);
        primeMpz(q);

        mpz_mul(n,p,q);    // n = pq

        // test
        // gmp_printf("\n\nC:%Zd\n",C);
        // gmp_printf("d:%Zd\n",d);
    
        // gmp_printf("r1:%Zd\n",r1);
        // gmp_printf("r2:%Zd\n",r2);
        // gmp_printf("r3:%Zd\n",r3);
        // gmp_printf("r4:%Zd\n",r4);

        // gmp_printf("p:%Zd\n",p);
        // gmp_printf("q:%Zd\n",q);
        // gmp_printf("n:%Zd\n\n\n",n);



        /*
         * Choose t1,t2,k in range 2~11 randomly.
         */
        t1 = random(2,11);
        t2 = random(2,11);
        k  = random(2,11);

        r  = random(2,11);
        t  = random(2,11);
        w  = random(2,11); 

        // cout<<"t1:"<<t1<<endl;             // Output t1,t2,k.
        // cout<<"t2:"<<t2<<endl;
        // cout<<"k :"<<k <<endl<<endl;



        /*
         * Preprocess
         */
        state = 1;                     // Initial the result state of mpz_invert()

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&beginT);

        // Calculate p1 = p^-1 mod q.
        state = mpz_invert(p1,p,q);    // If inverse exists return non-zero.
        if(state==0){
            cout<<"The inverse of p doesn't exist !"<<endl;
            exit(0);
        }

        // Calculate q1 = q^-1 mod p.
        state = mpz_invert(q1,q,p);
        if(state ==0){
            cout<<"The inverse of q doesn't exists !"<<endl;
            exit(0);
        }
    
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&endT);

        avgTimeA1 += timeCost(beginT,endT);



        /*
         * Problem transformation
         */

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&beginT);

        // Calculate dp,dq.
        mpz_sub_ui(p_1,p,1);
        mpz_mod(dp,d,p_1);

        mpz_sub_ui(q_1,q,1);
        mpz_mod(dq,d,q_1);

        // Calculate dp1,dp2.
        mpz_mul(tmp1,r1,p_1);
        mpz_add(dp1,dp,tmp1);

        mpz_mul(tmp1,r2,p_1);
        mpz_add_ui(tmp2,tmp1,k);
        mpz_mul_ui(tmp3,dp,t1);
        mpz_add(dp2,tmp2,tmp3);

        // Calculate dq1,dq2.
        mpz_mul(tmp1,r3,q_1);
        mpz_add(dq1,dq,tmp1);

        mpz_mul(tmp1,r4,q_1);
        mpz_add_ui(tmp2,tmp1,k);
        mpz_mul_ui(tmp3,dq,t2);
        mpz_add(dq2,tmp2,tmp3);

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&endT);
        avgTimeA2 += timeCost(beginT,endT);



        /*
         * Server Compute
         */

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&beginT);

        mpz_powm(R1 ,C,dp1,n);
        mpz_powm(RR1,C,dp2,n);
        mpz_powm(R2 ,C,dq1,n);
        mpz_powm(RR2,C,dq2,n);

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&endT);
        avgTimeA3 += timeCost(beginT,endT);
        


        /*
         * Verify
         * Already have used the powm() instead of pow().
         */

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&beginT);

        mpz_powm_ui(tmp1,R1,t1,p);
        mpz_powm_ui(tmp2,C,k,p);
        mpz_mul(tmp3,tmp1,tmp2);
        mpz_mod(tmp4,tmp3,p);
        mpz_mod(tmp5,RR1,p);
        cmp1 = mpz_cmp(tmp4,tmp5);    // Return "+" if tmp3>tmp4,"0" if ==,"-" if <

        mpz_powm_ui(tmp1,R2,t2,q);
        mpz_powm_ui(tmp2,C,k,q);
        mpz_mul(tmp3,tmp1,tmp2);
        mpz_mod(tmp4,tmp3,q);
        mpz_mod(tmp5,RR2,q);
        cmp2 = mpz_cmp(tmp4,tmp5);    // Return "+" if tmp3>tmp4,"0" if ==,"-" if <
        
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&endT);

        avgTimeA4 += timeCost(beginT,endT);

        if(cmp1!=0||cmp2!=0){
            cout<<"Plan A Verify failed!"<<cmp1<<","<<cmp2<<endl<<flush;
        }


        /*
         * Recovery
         */

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&beginT);

        mpz_mul(tmp1,p,p1);
        mpz_mul(tmp2,tmp1,R1);
        mpz_mul(tmp3,q,q1);
        mpz_mul(tmp4,tmp3,R2);
        mpz_add(tmp5,tmp2,tmp4);
        mpz_mod(M,tmp5,n);

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&endT);

        avgTimeA5 += timeCost(beginT,endT);



        /*
         * Direct compute by user.
         */
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&beginT);

        mpz_sub_ui(p_1,p,1);
        mpz_mod(dp,d,p_1);

        mpz_sub_ui(q_1,q,1);
        mpz_mod(dq,d,q_1);

        mpz_powm(Mp,C,dp,p);
        mpz_powm(Mq,C,dq,q);

        mpz_mul(tmp1,p,p1);
        mpz_mul(tmp2,tmp1,Mp);

        mpz_mul(tmp3,q,q1);
        mpz_mul(tmp4,tmp3,Mq);

        mpz_add(tmp5,tmp2,tmp4);

        mpz_mod(tmp6,tmp5,n);

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&endT);

        avgTimeA6 += timeCost(beginT,endT);



        /*
         * Another plan.
         * 
         */
        int cmpB;

        // T
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&beginT);

        mpz_mul(L,n,p);
        
        mpz_mul_ui(tmp1,n,r);
        mpz_add(CC,C,tmp1);

        //mpz_sub_ui(n_1,n,1);
        mpz_mul(n_1,p_1,q_1);
        mpz_mul(tmp1,k1,n_1);
        mpz_add(dd1,d,tmp1);

        mpz_mul_ui(tmp1,d,t);
        mpz_mul(tmp2,k2,n_1);
        mpz_add(tmp3,tmp1,tmp2);
        mpz_add_ui(dd2,tmp3,w);

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&endT);

        avgTimeB1 += timeCost(beginT,endT);

        // C
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&beginT);

        mpz_powm(Y1,CC,dd1,L);
        mpz_powm(Y2,CC,dd2,L);

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&endT);

        avgTimeB2 += timeCost(beginT,endT);

        // V
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&beginT);

        mpz_powm_ui(tmp1,C,w,n);
        mpz_powm_ui(tmp2,Y1,t,n);
        mpz_mul(tmp3,tmp1,tmp2);
        mpz_mod(tmp4,tmp3,n);

        mpz_mod(tmp5,Y2,n);

        cmpB = mpz_cmp(tmp4,tmp5);
        
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&endT);

        avgTimeB3 += timeCost(beginT,endT);

        if(cmpB!=0){
            cout<<"Plan B Verify failed!"<<cmpB<<endl<<flush;
        }



    }
    cout<<endl<<endl;;



    /*
     * Plan A
     */

    // Calculate the averge time of every phase.
    avgTimeA1 /= (double)round;
    avgTimeA2 /= (double)round;
    avgTimeA3 /= (double)round;
    avgTimeA4 /= (double)round;
    avgTimeA5 /= (double)round;
    avgTimeA6 /= (double)round;

    printf("Plan A\n");
    printf("    Preprocess             : %0.10lfs\n",avgTimeA1);
    printf("    Problem transformation : %0.10lfs\n",avgTimeA2);
    printf("    Compute                : %0.10lfs\n",avgTimeA3);
    printf("    Verify                 : %0.10lfs\n",avgTimeA4);
    printf("    Recovery               : %0.10lfs\n\n",avgTimeA5);

    // Calculate the user time cost and the server time cost.
    double userTimeA   = 0.0;
    double serverTimeA = 0.0;
    userTimeA   = avgTimeA2 + avgTimeA4 + avgTimeA5;
    serverTimeA = avgTimeA3;
    printf("    userTime               : %0.10lfs\n",userTimeA);
    printf("    DirectTime             : %0.10lfs\n",avgTimeA6);
    printf("    serverTime             : %0.10lfs\n\n\n",serverTimeA);



    /*
     * Plan B
     */
    avgTimeB1 /= (double)round;
    avgTimeB2 /= (double)round;
    avgTimeB3 /= (double)round;

    printf("Plan B\n");
    printf("    Problem transformation : %0.10lfs\n",avgTimeB1);
    printf("    Compute                : %0.10lfs\n",avgTimeB2);
    printf("    Verify                 : %0.10lfs\n\n",avgTimeB3);

    // Calculate the user time cost and the server time cost.
    double userTimeB   = 0.0;
    double serverTimeB = 0.0;
    userTimeB   = avgTimeB1 + avgTimeB3;
    serverTimeB = avgTimeB2;
    printf("    userTime               : %0.10lfs\n",userTimeB);
    printf("    serverTime             : %0.10lfs\n",serverTimeB);



    /*
     * Clear the memory of mpz.
     */
    mpz_clear(L);

    mpz_clear(k1);

    mpz_clear(d);
    mpz_clear(dd1);
    mpz_clear(dd2);

    mpz_clear(M);
    mpz_clear(C);
    mpz_clear(CC);

    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(n);

    mpz_clear(p1);
    mpz_clear(q1);
    
    mpz_clear(dp);
    mpz_clear(dq);

    mpz_clear(Mp);
    mpz_clear(Mq);

    mpz_clear(Y1);
    mpz_clear(Y2);

    mpz_clear(p_1);
    mpz_clear(q_1);
    mpz_clear(n_1);

    //mpz_clear(r);
    mpz_clear(r1);
    mpz_clear(r2);
    mpz_clear(r3);
    mpz_clear(r4);

    mpz_clear(dp1);
    mpz_clear(dp2);
    mpz_clear(dq1);
    mpz_clear(dq2);

    mpz_clear(R1);
    mpz_clear(RR1);
    mpz_clear(R2);
    mpz_clear(RR2);

    mpz_clear(tmp1);
    mpz_clear(tmp2);
    mpz_clear(tmp3);
    mpz_clear(tmp4);
    mpz_clear(tmp5);
    mpz_clear(tmp6);

    return 0;
}
