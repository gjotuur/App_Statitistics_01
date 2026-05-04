#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <stdlib.h>
#include <openssl/rand.h>

/*To compile: gcc main.c -o tst1.exe -IC:/msys64/mingw64/include -LC:/msys64/mingw64/lib -lcrypto -lssl*/

double get_rand_01(){
    uint64_t value;

    RAND_bytes((unsigned char*)&value, sizeof(value));
    return (double)value / (double)UINT64_MAX;
}

typedef struct{
    uint64_t size;
    double* samples;
}stat_sample_t;

void Init_Sample(stat_sample_t* s, uint64_t size){
    s->size = size;
    s->samples = calloc(s->size,sizeof(double));
}

void Fill_Sample_01(stat_sample_t* s){
    for(uint64_t i = 0; i < s->size; i++){
        s->samples[i] = get_rand_01();
    }
}

int main(){
    stat_sample_t* sample1 = malloc(sizeof(stat_sample_t));
    Init_Sample(sample1, (uint64_t)100);
    Fill_Sample_01(sample1);

    for(uint64_t i = 0; i < sample1->size; i++){
        printf("omega_%d = %.2lf, ", i+1, sample1->samples[i]);
        ((i+1)%10 == 0) ? printf("\n") : 0;
    }
    double gay_website = get_rand_01();
    printf("\nGay website is located %.2lf blocks down\n", gay_website);

    return 0;
}
