#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <stdlib.h>
#include <openssl/rand.h>
#include <string.h>
#include <stdbool.h>
#define PI 3.14159265358979323946

/*To compile: gcc main.c -o tst1.exe -IC:/msys64/mingw64/include -LC:/msys64/mingw64/lib -lcrypto -lssl*/


//Basic struct for samples (w/o mean and dispersion)
typedef struct{
    uint64_t size;
    double* samples;
}stat_sample_t;

//Block of funcs used to work with samples
void Init_Sample(stat_sample_t* s, uint64_t size);                                  //Sample struct initialization with size, use malloc before calling
stat_sample_t* Copy_Sample(stat_sample_t* s1);                                      //Copying samples, memory allocation integrated, return stat_sample_t* value, example of use: stat_sample_t* s2 = Copy_Sample(s1);
void Fill_Sample_with_random(stat_sample_t* s);                                     //Filling the sample with pseudo-random numbers from [0;1] interval
void Normalize_Sample(stat_sample_t* s);                                            //Setting sample elements to fit with normal distribution N(0,1)
double Sample_Mean(stat_sample_t* s);                                               //Sample mean evaluation
double Sample_Variance(stat_sample_t* s, bool bias);                                //Sample variance evaluation, 2nd param set 1 for unbiased and 0 for biased
void Sample_Analyze(stat_sample_t* s, double* mean, double* var,                    //Single function to evaluate all at once, params is pointers to store values
                    bool bias, double* st_dev);


int main(){
    stat_sample_t* sample1 = malloc(sizeof(stat_sample_t));
    Init_Sample(sample1, (uint64_t)10000);
    Fill_Sample_with_random(sample1);

    for(uint64_t i = 0; i < sample1->size; i++){
        printf("omega_%d = %.2lf, ", i+1, sample1->samples[i]);
        ((i+1)%10 == 0) ? printf("\n") : 0;
    }

    stat_sample_t * sample2 = Copy_Sample(sample1);

    Normalize_Sample(sample2);

    for(uint64_t i = 0; i < sample2->size; i++){
        printf("ksi_%d = %.2lf, ", i + 1, sample2->samples[i]);
        ((i+1)%10 == 0) ? printf("\n") : 0;
    }

    double mean = Sample_Mean(sample1);
    double var = Sample_Variance(sample1, true);

    double mean_p, var_p, st_dev_p;
    Sample_Analyze(sample1, &mean_p, &var_p, true, &st_dev_p);

    printf("\nResults of single functions:\nMean = %.2lf   Variance = %.2lf    St_Dev = %.2lf", mean, var, sqrt(var));
    printf("\nResults of multi function:\nMean = %.2lf   Variance = %.2lf    St_Dev = %.2lf", mean_p, var_p, st_dev_p);

    return 0;
}

//Utility pseudo-random number generation function using OpenSSL high-quality randomness
double get_rand_01(){
    uint64_t value;

    RAND_bytes((unsigned char*)&value, sizeof(value));
    return (double)value / (double)UINT64_MAX;
}


void Init_Sample(stat_sample_t* s, uint64_t size){
    s->size = size;
    s->samples = calloc(s->size,sizeof(double));
}

stat_sample_t* Copy_Sample(stat_sample_t* s1){

    if(s1 == NULL) return NULL;

    stat_sample_t* s2 = malloc(sizeof(stat_sample_t));
    if(s2 == NULL) return NULL;
    
    Init_Sample(s2, s1->size);
    if(s2->samples == NULL){
        free(s2);
        return NULL;
    }

    memcpy(s2 -> samples, s1 -> samples, s1->size * sizeof(double));

    return s2;
}


void Fill_Sample_with_random(stat_sample_t* s){
    for(uint64_t i = 0; i < s->size; i++){
        s->samples[i] = get_rand_01();
    }
}


void Normalize_Sample(stat_sample_t* s){
    double ksi_1 = 0;
    double ksi_2 = 0;

    for(uint64_t i = 0; i < s->size; i += 2){

        if(i+1 <= s->size){
            ksi_1 = sqrt(-2*log(s->samples[i]))*sin(2*PI*s->samples[i+1]);
            ksi_2 = sqrt(-2*log(s->samples[i]))*cos(2*PI*s->samples[i+1]);

            s->samples[i] = ksi_1;
            s->samples[i+1] = ksi_2;

        } else {
            double omega_plus = get_rand_01();                                                              //if sample size is even, just add 1 more random value to evaluate last element

            ksi_1 = sqrt(-2*log(s->samples[i]))*sin(2*PI*omega_plus);

            s->samples[i] = ksi_1;
        }
    }
}

double Sample_Mean(stat_sample_t* s){

    if (s == NULL || s->size == 0) return 0;

    double s_mean = 0;

    for(uint64_t i = 0; i < s->size; i++){
        s_mean += s->samples[i];
    }
    return (s_mean / (double)s->size);
}

double Sample_Variance(stat_sample_t* s, bool bias){
    if(s == NULL || s-> size < 2) return 0;

    double s_var = 0;
    double bias_q = (bias) ? 1 : 0;
    double mean = Sample_Mean(s);

    for(uint64_t i = 0; i < s->size; i++){
        double diff = s->samples[i] - mean;
        s_var += diff*diff;
    }

    return (s_var / (double)(s->size - bias_q));
}

void Sample_Analyze(stat_sample_t* s, double* mean, double* var, bool bias, double* st_dev){

    if(s == NULL) return;
    if(s->size < 2 || s->size == 0) return;

    double s_mean = 0;
    double s_var = 0;
    uint64_t q_bias = (bias) ? 1 : 0;

    for(uint64_t i = 0; i < s->size ; i++){
        s_mean += s->samples[i];
    }

    s_mean /= (double)s->size;

    for(uint64_t i = 0; i < s->size; i++){
        double diff = s->samples[i] - s_mean;
        s_var += diff * diff;
    }

    s_var /= (double)(s->size - q_bias);
    double s_st_dev = sqrt(s_var);

    *mean = s_mean;
    *var = s_var;
    *st_dev = s_st_dev;
}