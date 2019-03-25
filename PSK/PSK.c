#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>
#include<gsl/gsl_randist.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#define PI M_PI
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
//#define NoE 100

int N;
void reverse(char str[], int length,char s[]) 
{ 
    int start = 0; 
    int end = length -1; 
    char a, b, temp;
    int l = end;
    for(int i = start; i <= (end-start+1)/2;i++){

        a = str[i];
        b = str[l];
        temp = b;
        b = a;
        a = temp;
        s[i] = a;
        s[l] = b;
        l--;
    }
    s[end+1]='\0';
} 

int bintogray(int bin)
{
    int a, b, result = 0, i = 0;
    while (bin!=0)
    {
        a = bin % 10;
        bin /= 10;
        b = bin % 10;
        if ((a && !b) || (!a && b))
            result += pow(10, i);
        i++;
    }
    return result;
}


int generateCodes(int n, char gray[][25]){
	char A[n][25];
    int num, i, j;
    int k = log(n)/log(2), len = k;
    if (n == 2){
    	gray[0][0] = '0';gray[0][1] = '\0';
    	gray[1][0] = '1';gray[1][1] = '\0';
    	return k;
    }
    for(i = 0; i < n; i++){
        char str[k+1];
        char s[k+1];
        num = i;
        for(int j=len-1;j>=0;j--){
            A[i][j] = num%2+'0';
            num /= 2;
        }
        A[i][len]='\0';
        num = bintogray(atoi(A[i]));
        int l =0;
        if (num == 0) 
        	str[l++] = '0';
        while (num != 0) 
        { 
            int rem = num % 10; 
            str[l++] = (rem > 9) ? (rem-10) + 'a' : rem + '0'; 
            num /= 10; 
        } 
        if(strlen(str) < k){            
            for(j = strlen(str); j < k; j++)
                str[j] = '0';            
        }
        str[j] = '\0'; 
        reverse(str,k,s);
        strcpy(gray[i],s);
    }
    return k;
}

void generatePSK(float cons[][2], float angle[], int N, float R){
	int b;
	for(b = 0; b < N; b++)
    {
        angle[b] = (b*2.0*PI)/N;
        cons[b][0] = R * cos(angle[b]);
        cons[b][1] = R * sin(angle[b]);
    }
}

float raisedCosine(float *rc){
    float T=1.0,a=0.25;
    float t=0.01,Eg,temp;

    double data[2*700];
    FILE *fft;
    fft = fopen("fft.txt","w");  
    
    for(int t1=0;t1<700;t1++){
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
        t += 0.01;
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
    if(N == 32)
    	return Eg/2;
    return Eg*1.5;
}

int main()
{
    gsl_rng_env_setup();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);

    printf("Enter the value of N(4,8,16,32,64)= ");
    scanf("%d",&N);

    FILE *qc;
    qc = fopen("con.txt","w");

    int NoE = (N == 4 || N == 2) ? 10 : 100;
	
    char gray[N][25];
    int k = generateCodes(N, gray);
    float angle[N];
    float cons[N][2];
    float R = 1;
    generatePSK(cons, angle, N, R);
    for(int b = 0; b < N; b++)
    {
        printf("(%f, %f)	%s\n",cons[b][0], cons[b][1], gray[b]);
        fprintf(qc, "%f %f\n",cons[b][0],cons[b][1]);
    }
    fclose(qc);
    int i, j, count, iter, n;
    
    float *rc = (float *)malloc(1000*sizeof(float));
    float Eg = raisedCosine(rc);
    //printf("E = %f\n",Eg);
    
	FILE *fp = fopen("PSK.txt","a");
    float sm[2001], dummy_f[2001][2], N0, std_noise;
    int f = 20;
    for(int snr = 1; snr <= 10; snr++)
    {
        count = 0;
        N0=(2*Eg)/snr;
        std_noise = sqrt(N0/2);
        char bits[k+1];
        float psk[2], thetha; 
        iter = 0;
        while(iter < 3000000)
        //for(iter=0;;iter++)
        {
            for(i = 0; i < k; i++)
            {
                bits[i] = gsl_rng_uniform_int(r,2) + '0';
            }
            bits[i] = '\0';
            for(j = 0;j < N;j++)
            {
                if(strcmp(bits,gray[j]) == 0)
                {
                    psk[0] = cons[j][0];
                    psk[1] = cons[j][1];
                    thetha = angle[j];
                    break;
                }
            }
            float gn=gsl_ran_gaussian(r,std_noise);
            float t = 0.01;
            for(int t1=0;t1<700;t1++)
            {
                //float gn=gsl_ran_gaussian(r,std_noise);
                dummy_f[t1][0] = sqrt(2/Eg) * rc[t1] * cos(2*PI*f*t);
                dummy_f[t1][1] = -1 * sqrt(2/Eg) * rc[t1] * sin(2*PI*f*t);
            
                sm[t1]=(sqrt(Eg/2)*cos(thetha)*dummy_f[t1][0]) + (sqrt(Eg/2)*sin(thetha)*dummy_f[t1][1]) + gn;
                t += 0.01;
            }
            float x0 = 0,xn = 1,h = 0.01 , so1 = 0, so2 = 0, ans_x, ans_y;
            n=(xn-x0)/h;
            if(n%2==1)
            	n++;
            h=(xn-x0)/n;
            for(i=0; i<n; i++)
            {
                so1 += sm[i]*dummy_f[i][0];
                so2 += sm[i]*dummy_f[i][1];
            }
            int temp = 1;
            if(N == 16)
                temp = 15;
            ans_x= so1/temp;
            ans_y=so2/temp;
            //printf("\nfinal integration is :  %f %f\n\n",ans_x,ans_y);
            float d, min = FLT_MAX;
            float est_x, est_y;
            int min_at;
            for (i = 0;i < N; i++){
                d = sqrt(pow(ans_x - cons[i][0], 2) + pow(ans_y - cons[i][1], 2));
                if (min > d){
                    min_at = i;
                    min = d;
                }
            }
            est_x = cons[min_at][0];est_y = cons[min_at][1];
            
            if(psk[0]!=est_x || psk[1]!=est_y)
                count++;
            /*else
                printf("(%f, %f) ------------ (%f, %f)\n",psk[0],psk[1],est_x,est_y);*/
            if(count == NoE)
                break;
            //if(iter > 300000)
            //	break;
        iter ++;
        }
        fprintf(fp, "%d %0.14f\n",snr,(float)count/iter);
        printf("**********count = %d in iterations = %d at SNR = %d dB\n",count,iter,snr);
        printf("**********probability of error is %0.14f\n", (float)count/iter);
    }
    fprintf(fp,"\n\n");
    fclose(fp);
}
