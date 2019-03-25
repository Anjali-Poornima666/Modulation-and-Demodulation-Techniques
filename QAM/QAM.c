#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>
#include<gsl/gsl_randist.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
#define PI M_PI
#define NoE 100


void generateCodes(char A[][7], int n, int len){
    for(int i = 0; i < n; i++){
        int num = i;
        for(int j = len-1; j >= 0; j--){
            A[i][j]=num%2+'0';
            num=num/2;
        }
        A[i][len]='\0';
    }
}

void generateQAM(int cons[][2], int qam_n){
	int x = (sqrt(qam_n)-floor(sqrt(qam_n)) == 0.0) ? (int)sqrt(qam_n) : (int)sqrt(qam_n*2) ;
	int y = qam_n/x, c = 0;
	int ref1 = -1*(x-1), ref2 = -1*(y-1), flag1 = 0, flag2 = 0;
	for(int i=0;i<qam_n;i++){
		if(c == y){
			ref1 = (flag1 == 0) ? ref1 + 2 : ref1 - 2;
			c=0;
			if(ref1==1 && flag1==0){
				ref1 = x-1;
				flag1 = 1;
			}
		}
		cons[i][0] = ref1;
		cons[i][1] = ref2;
		ref2 = (flag2 == 0) ? ref2 + 2 : ref2 - 2;
		if(ref2 == -1 && flag2 == 1){
			flag2=0;
			ref2=-1*(y-1);
		}
		else if((ref2 == 1 || ref2 == 2) && flag2 == 0){
			flag2=1;
			ref2=y-1;
		}
		c++;
	}
}

float raisedCosine(float *rc){
    float T=1.0,a=0.25;
    float t=0.01,Eg,temp;

    double data[2*128];
    FILE *fft;
    fft = fopen("fft.txt","w");  
    
    for(int t1=0;t1<128;t1++){
        if (t1 == (T/(2.0*a))*100 - 1){
    		*(rc+t1) = (1/(4.0*T))*sin(PI/(2.0*a))*2.0*a;
    	}
    	else if (t == 0)
    		*(rc+t1) = 1.0;
    	else{
		    *(rc+t1) = (1.0/(PI*t))*sin(PI*t/T)*cos(PI*a*t/T);
		    temp = 1 - pow((2.0*a*t/T), 2);
		    *(rc+t1) = *(rc+t1)/temp;
		}
		Eg += pow(rc[t1],2);
		//printf("%f %f\n", *(rc + t1), t);
		REAL(data,t1) = rc[t1];
        IMAG(data,t1) = 0.0;
        t=t+0.01;
    }

    gsl_fft_complex_radix2_forward (data, 1, 128);
    t = 0.01;
    for (int i = 0; i < 128; i++)
    {
        fprintf(fft,"%f %e\n",t,REAL(data,i));
        //printf ("%d %e %e\n", i,REAL(data,i),IMAG(data,i));
        t += 0.1;
    }
    fclose(fft);
    //printf("asd\n");
    return Eg/2;
}

int main()
{
    gsl_rng_env_setup();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);
    
    int N, i, j, count, iter;
    printf("Enter the value of N(4,8,16,32,64)= ");
    scanf("%d", &N);
    int k = log10(N)/log10(2);
    char A[N][7];
    generateCodes(A, N, k);
    int cons[N][2];
    generateQAM(cons, N);

    FILE *qc;
    qc = fopen("con.txt","w");

    for(i = 0;i < N;i++){
        printf("%s : %d+%d*j\n",A[i],cons[i][0],cons[i][1]);
        fprintf(qc, "%d %d\n",cons[i][0],cons[i][1]);
    }
    fclose(qc);
    //srand(time(NULL));

    float *rc = (float *)malloc(1000*sizeof(float));
    float Eg = raisedCosine(rc);
    //printf("E = %f\n",Eg);
    printf("*********************** %d QAM *******************************\n",N);

    float sm[2001], dummy_f[2001][2], N0, std_noise;
    int f = 20, qam[2];
    char bits[k+1];

    FILE *fp = fopen("QAM.txt","w");

    for(int snr=1;snr<=10;snr++)
    {
        count=0;
        N0=(2*Eg)/snr;
        std_noise = sqrt(N0/2);
        iter =0;
        while(iter < 300000)        
        //for(iter=0;;iter++)
        {
            for(i=0;i<k;i++)
            {
                bits[i] = gsl_rng_uniform_int(r,2) + '0';
            }
            bits[i] = '\0';
            for(j = 0;j < N;j++)
            {
                if(strcmp(bits, A[j]) == 0)
                {
                    qam[0] = cons[j][0];
                    qam[1] = cons[j][1];
                    break;
                }
            }
            //printf("%s %d  %d\n",bits, qam[0], qam[1]);
            //float gn = gsl_ran_gaussian(r,std_noise);
            float t = 0.01;
            for(int t1 = 0; t1 < 128; t1++)
            {
            	float gn = gsl_ran_gaussian(r,std_noise);
                dummy_f[t1][0] = sqrt(2/Eg) * rc[t1] * cos(2*PI*f*t);
                dummy_f[t1][1] = -1 * sqrt(2/Eg) * rc[t1] * sin(2*PI*f*t);
            
                sm[t1] = (sqrt(Eg/2)*qam[0]*dummy_f[t1][0]) + (sqrt(Eg/2)*qam[1]*dummy_f[t1][1]) + gn;
                t += 0.01;
            }
            int n;
            float x0 = 0, xn = 1, h = 0.01, so1 = 0, so2 = 0, ans_x, ans_y;
            n=(xn-x0)/h;
            if(n%2==1)
                n++;
            h=(xn-x0)/n;
            for(i=0; i<n; i++)
            {
                so1 += (sm[i]*dummy_f[i][0]);
                so2 += (sm[i]*dummy_f[i][1]);
            }
            ans_x = so1/6.2;
            ans_y = so2/6.2;
            //printf("\nfinal integration is :  %f %f\n\n",ans_x,ans_y);
            float d, min = FLT_MAX;
            int est_x, est_y, min_at;
            for (i = 0;i < N; i++){
                d = sqrt(pow(ans_x - cons[i][0], 2) + pow(ans_y - cons[i][1], 2));
                if (min > d){
                    min_at = i;
                    min = d;
                }
            }
            est_x = cons[min_at][0];est_y = cons[min_at][1];
            //printf("(%d, %d)		(%d, %d)\n",qam[0],qam[1],est_x,est_y);
            
            if(qam[0]!=est_x || qam[1]!=est_y)
                count++;
            if(count == NoE)
                break;
            //if(iter >= 300000)
            	//break;
        iter++;
        }
        fprintf(fp, "%d %0.14f\n",snr,(float)count/iter);
        printf("**********count = %d in iterations = %d at SNR = %d dB\n",NoE,iter,snr);
        printf("**********probability of error is %.14f\n", (float)count/iter);
    }
    fclose(fp);
}
