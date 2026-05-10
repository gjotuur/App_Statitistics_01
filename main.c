#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <stdlib.h>
#include <openssl/rand.h>
#include <string.h>
#define PI 3.14159265358979323946

/*To compile: gcc main.c -o tst1.exe -IC:/msys64/mingw64/include -LC:/msys64/mingw64/lib -lcrypto -lssl*/


//Pseudo-random numbers generation using OpenSSL high-quality randomness
double get_rand_01(){
    uint64_t value;

    RAND_bytes((unsigned char*)&value, sizeof(value));
    return (double)value / (double)UINT64_MAX;
}

//Basic struct for samples (w/o mean and dispersion)
typedef struct{
    uint64_t size;
    double* samples;
}stat_sample_t;


//Sample struct initialization func, use malloc for stat_sample before calling
void Init_Sample(stat_sample_t* s, uint64_t size){
    s->size = size;
    s->samples = calloc(s->size,sizeof(double));
}

//Make the deep copy of sample using memcpy func
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


//Use pseudo-random generator to fill the sample with values
void Fill_Sample_with_random(stat_sample_t* s){
    for(uint64_t i = 0; i < s->size; i++){
        s->samples[i] = get_rand_01();
    }
}


//Set a values inside sample to conform to a normal distribution
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

int main(){
    stat_sample_t* sample1 = malloc(sizeof(stat_sample_t));
    Init_Sample(sample1, (uint64_t)100);
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

    return 0;
}
