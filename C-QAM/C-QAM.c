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

void generateQAM(int NoC, float cons[][2], char gray[][25]){

    FILE *qc;
    qc = fopen("con.txt","w");

	float angle[2] = {0, PI/4};
	int radius, count = 1;
    for(int b = 0; b < N; b++)
    {
        for(int rem1 = 0; rem1 < NoC; rem1++){
            if(b%NoC == rem1)
           	{
                radius = rem1+1;
                break;
            }
        }
        int a = radius%2;
        cons[b][0] = radius * cos(angle[a]);
        cons[b][1] = radius * sin(angle[a]);
        printf("%f    (%f, %f)   %s  \n",angle[a],cons[b][0],cons[b][1],gray[b]);
        fprintf(qc, "%f %f\n",cons[b][0],cons[b][1]);
        if(count == NoC){
        	count = 0;
        	angle[0] += PI/2;
        	angle[1] += PI/2;
        }
        count++;
    }
    fclose(qc);
}

float raisedCosine(float *rc){
    float T=1.0,a=0.25;
    float t=0.01,Eg,temp;
    
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
        t += 0.01;
    }
    return Eg/2;
}

int main()
{
	gsl_rng_env_setup();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);
	printf("Enter the value of N(8,16)= ");
    scanf("%d",&N);
    
    char gray[N][25];
    int NoC;
    if(N == 2)
    	NoC = 1;
   	else
   		NoC = N/4;
    int k = generateCodes(N,gray);
    float cons[N][2];
    generateQAM(NoC, cons, gray);
    int i, j, iter=0, count;
    
    srand(time(NULL));
    
    float *rc = (float*)malloc(1000*sizeof(float));
    float Eg = raisedCosine(rc);
    //printf("E = %f\n",Eg);
	FILE *fp = fopen("C-QAM.txt", "a");
    float sm[2001],dummy_f[2001][2],N0;
    int f = 20;
    float std_noise,temp;
    for(int snr=1;snr<=10;snr++)
    {
        count = 0;
        N0 = (2*Eg)/snr;
        std_noise = sqrt(N0/2);
        char bits[k+1];
        float qam[2];
        iter = 0;
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
                if(strcmp(bits,gray[j]) == 0)
                {
                    qam[0] = cons[j][0];
                    qam[1] = cons[j][1];
                    break;
                }
            }
            
            float t=0.01;
            for(int t1=0;t1<700;t1++)
            {
            	float gn=gsl_ran_gaussian(r,std_noise);
                dummy_f[t1][0]=sqrt(2/Eg) * rc[t1] * cos(2*PI*f*t);
                dummy_f[t1][1]=-1*sqrt(2/Eg) * rc[t1] * sin(2*PI*f*t);
                
                sm[t1]=(sqrt(Eg/2)*qam[0]*dummy_f[t1][0]) + (sqrt(Eg/2)*qam[1]*dummy_f[t1][1]) + gn;
                t=t+0.01;
            }
            int n;
            float x0 = 0,xn = 1,h = 0.01 ,so1 = 0, so2 = 0,ans_x,ans_y;
            n=(xn-x0)/h;
            if(n%2==1)
            	n++;
            h=(xn-x0)/n;
            for(i=0; i<n; i++)
            {
                so1 += sm[i]*dummy_f[i][0];
                so2 += sm[i]*dummy_f[i][1];
            }
            float temp = 6.2;
            ans_x= so1/temp;
            ans_y=so2/temp;
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
            
            if(qam[0]!=est_x || qam[1]!=est_y)
                count++;
            if(count == NoE)
                break;
			//if(iter>300000)
            //	break;
            iter ++;
        }
        fprintf(fp, "%d %0.14f\n",snr,(float)count/iter);
        printf("**********count = %d in iterations = %d at SNR = %d dB\n",count,iter,snr);
        printf("**********probability of error is %0.14f\n", (float)count/iter);
    }
  	fprintf(fp, "\n\n");
  	fclose(fp);
}
