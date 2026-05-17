#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <stdlib.h>
#include <openssl/rand.h>
#include <string.h>
#include <stdbool.h>
#include <stdarg.h>
#define PI 3.14159265358979323946

/*To compile: gcc main.c -o tst1.exe -IC:/msys64/mingw64/include -LC:/msys64/mingw64/lib -lcrypto -lssl*/

//type to store lookup tables for distributions
typedef struct{
    double alpha;
    uint64_t n;
    double value;
}Distribution_t;

//Struct to store names of used distributions
typedef enum {
    Student = 0,
    Normal = 1,
    Khi_low = 2,
    Khi_high = 3
}Distribution_name_t;

typedef enum{
    Weibull_Ksi = 0,
    Weibull_EtaA = 1,
    Weibull_EtaB = 2
} Weibull_Type_t;

//Basic struct for samples (w/o mean and dispersion)
typedef struct{
    uint64_t size;
    double* samples;
}stat_sample_t;

//Confidence intervals always include 2 double values
typedef struct{
    double upper_limit;
    double lower_limit;
}Confidence_Interval;

//Utility structure to store results
typedef struct{
    uint64_t size;
    Confidence_Interval method1;
    Confidence_Interval method2;
    Confidence_Interval method3;
    double mean;
    double var;
    double st_dev;  
}Sample_Test_Result;

//Utility structure to store results of part 2
typedef struct{
    uint64_t final_N;
    double Q;
    double variance;
    Confidence_Interval cnf;
}convergence_result_t;


//Block of funcs used to work with numbers (generation, transformation)
double get_rand_01();                                                               //Returns pseudo-random number, uses high-quality randomness from OpenSSL
double Set_Number_To_Weibull(double n, Weibull_Type_t type);                        //Generate random [0,1] or use get_rand_01, function transform value to Weibull distribution


//Block of funcs used to work with samples
void Init_Sample(stat_sample_t* s, uint64_t size);                                  //Sample struct initialization with size, use malloc before calling
void Clear_Samples(stat_sample_t* first, ...);                                      //Clear all the samples used, list always must end with NULL-terminator
stat_sample_t* Copy_Sample(stat_sample_t* s1);                                      //Copying samples, memory allocation integrated, return stat_sample_t* value, example of use: stat_sample_t* s2 = Copy_Sample(s1);
void Fill_Sample_with_random(stat_sample_t* s);                                     //Filling the sample with pseudo-random numbers from [0;1] interval
void Normalize_Sample(stat_sample_t* s);                                            //Setting sample elements to fit with normal distribution N(0,1)
double Sample_Mean(stat_sample_t* s);                                               //Sample mean evaluation
double Sample_Variance(stat_sample_t* s, bool bias);                                //Sample variance evaluation, 2nd param set 1 for unbiased and 0 for biased
void Sample_Analyze(stat_sample_t* s, double* mean, double* var,                    //Single function to evaluate all at once, params is pointers to store values
                    bool bias, double* st_dev);
double Distribution_Value(uint64_t n, double gamma, Distribution_name_t d);         //Single function to find value of a distribution with a given parameters of gamma and n, uses lookup tables, input for distribution name may be {Normal; Student; Khi_...}
void Confidence_Interval_Var1(stat_sample_t* sample, double gamma,                  //Confidence interval for mean when variance is unknown and distribution is N(0,1)
                                double* lower_limit, double *upper_limit);
void Confidence_Interval_Var2(stat_sample_t* sample, double gamma,                  //Confidence interval for variance for sample with distribution type of N(0,1)
                                double* lower_limit, double* upper_limit);
void Confidence_Interval_Var3(stat_sample_t* sample, double gamma,                  //Confidence interval in case if the distribution is unknown
                                double* lower_limit, double* upper_limit);
void Set_Sample_To_Weibull(stat_sample_t* s, Weibull_Type_t type);                  //To use first generate sample with random values on [0,1] interval, copy if needed, func will refill existing sample with Weibull distribution values

//Convergence test
convergence_result_t Convergence_Test_V1(uint64_t sample_size, Weibull_Type_t variation);               //Monte-Carlo method
convergence_result_t Convergence_Test_V2(uint64_t sample_size, Weibull_Type_t variation);               //Conditional mean method




int main(){

    /*uint64_t sample_sizes[] = {10, 100, 1000, 10000, 100000, 1000000};
    int sample_quantity = 6;

    Sample_Test_Result* Sample_Results = calloc(sample_quantity, sizeof(Sample_Test_Result));

    printf("Task #1: Confidence intervals");

    double gamma = 0.01;

    for(int i = 0; i < sample_quantity; i++){
        stat_sample_t* sample = malloc(sizeof(stat_sample_t));
        Init_Sample(sample, sample_sizes[i]);

        Fill_Sample_with_random(sample);
        Normalize_Sample(sample);

        printf("Test %d, sample size n = %7llu, gamma = 0.01", i, sample->size);

        Sample_Results[i].size = sample_sizes[i];
        
        Sample_Analyze(sample, &Sample_Results[i].mean, &Sample_Results[i].var, true, &Sample_Results[i].st_dev);

        Confidence_Interval_Var1(sample, gamma, &Sample_Results[i].method1.lower_limit, &Sample_Results[i].method1.upper_limit);
        Confidence_Interval_Var3(sample, gamma, &Sample_Results[i].method2.lower_limit, &Sample_Results[i].method2.upper_limit);
        Confidence_Interval_Var2(sample, gamma, &Sample_Results[i].method3.lower_limit, &Sample_Results[i].method3.upper_limit);

        printf(" ......done\n");

        Clear_Samples(sample, NULL);
    }


    printf("\nResults are");
    for(int i = 0; i < sample_quantity; i++){
        printf("\nSample %2d: size = %7d, interval M1 [%.6lf, %.6lf], interval M2 [%.6lf, %.6lf], Interval M3 [%.6lf, %.6lf]", 
            i, Sample_Results[i].size, Sample_Results[i].method1.lower_limit, Sample_Results[i].method1.upper_limit,
            Sample_Results[i].method2.lower_limit, Sample_Results[i].method2.upper_limit,
            Sample_Results[i].method3.lower_limit, Sample_Results[i].method3.upper_limit);

        printf("\nReal mean = %.6lf, real variance = %.6lf, ranges is [%.6lf, %.6lf, %.6lf]", Sample_Results[i].mean, Sample_Results[i].var,
        (-1)*Sample_Results[i].method1.lower_limit + Sample_Results[i].method1.upper_limit,
        (-1)*Sample_Results[i].method2.lower_limit + Sample_Results[i].method2.upper_limit,
        (-1)*Sample_Results[i].method3.lower_limit + Sample_Results[i].method3.upper_limit);
    }*/

    stat_sample_t* sample1 = malloc(sizeof(stat_sample_t));

    Init_Sample(sample1,100);
    Fill_Sample_with_random(sample1);

    stat_sample_t* sample2 = Copy_Sample(sample1);
    Set_Sample_To_Weibull(sample2, Weibull_EtaB);

    for(uint64_t i = 0; i < sample1->size; i++){
        printf("   s1[%3llu] = %2.6lf | s2[%3llu] = %2.6lf", i+1, sample1->samples[i], i+1, sample2->samples[i]);
        //double s_i = sample1->samples[i];
        //double s_2i = Set_Number_To_Weibull(s_i, Weibull_Ksi);
        //printf("    s1[%3llu] = %2.6lf | s2[%3llu] = %2.6lf", i+1, s_i, i+1, s_2i);
        (i+1) % 5 == 0 ? printf("\n") : 0;
    }


    Clear_Samples(sample1, sample2, NULL);

    convergence_result_t try_1 = Convergence_Test_V1(200, Weibull_EtaA);
    printf("\nResults are \nN = %llu\nVariance = %.4lf\nQ = %lf\nInterval [%.4lf,%.4lf], range is %.4lf", try_1.final_N, try_1.variance, try_1.Q, try_1.cnf.lower_limit, try_1.cnf.upper_limit, try_1.cnf.upper_limit - try_1.cnf.lower_limit);

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

void Set_Sample_To_Weibull(stat_sample_t* s, Weibull_Type_t type){

    if(s == NULL || s->size == 0) return;

    double omega;
    switch(type){
        case Weibull_Ksi:
            for(uint64_t i = 0; i < s->size; i++){
                omega = sqrt(-log((s->samples[i] > 0 ? s->samples[i] : 1e-15)));
                s->samples[i] = omega;
            }
            break;
        case Weibull_EtaA:
            for(uint64_t i = 0; i < s->size; i++){
                omega = (double)1 / sqrt((s->samples[i] > 0 ? s->samples[i] : 1e-15)) - 1;
                s->samples[i] = omega;
            }
            break;
        case Weibull_EtaB:
            for(uint64_t i = 0; i < s->size; i++){
                omega = -log(s->samples[i] > 0 ? s->samples[i] : 1e-15);
                s->samples[i] = omega;
            }
    }
}

double Set_Number_To_Weibull(double n, Weibull_Type_t type){
    double converted;
    if(type == Weibull_Ksi) converted = sqrt(-log((n > 0 ? n : 1e-15)));
    if(type == Weibull_EtaA) converted = (double)1 / sqrt((n > 0 ? n : 1e-15)) - 1;
    if(type == Weibull_EtaB) converted = -log(n > 0 ? n : 1e-15);

    return converted;
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

void Clear_Samples(stat_sample_t* first, ...){
    if(first == NULL) return;

    va_list args;
    va_start(args, first);

    stat_sample_t* stream = first;

    while(stream != NULL){
        if(stream->samples != NULL){
            free(stream->samples);
            stream->samples = NULL;
        }
        free(stream);

        stream = va_arg(args, stat_sample_t*);
    }

    va_end(args);
}

double Distribution_Value(uint64_t n, double gamma, Distribution_name_t d){

    double alpha = gamma;
    if(d == Normal){
        if(alpha == 0.05) return 1.9600;
        if(alpha == 0.01) return 2.5758;
        if(alpha == 0.001) return 3.2905;
        if(alpha == 0.1) return 1.6449;
        return 1.9600;
    }

    static const Distribution_t Student_tab[] = {
        /*n == 10*/     {0.1,      10, 1.833}, {0.05,      10, 2.262}, {0.01,      10, 3.250}, {0.001,      10, 4.781},
        /*n == 100*/    {0.1,     100, 1.660}, {0.05,     100, 1.984}, {0.01,     100, 2.626}, {0.001,     100, 3.390},
        /*n == 1000*/   {0.1,    1000, 1.646}, {0.05,    1000, 1.962}, {0.01,    1000, 2.581}, {0.001,    1000, 3.300},
        /*n == 10000*/  {0.1,   10000, 1.645}, {0.05,   10000, 1.960}, {0.01,   10000, 2.576}, {0.001,   10000, 3.291},
        /*n == 100000*/ {0.1,  100000, 1.645}, {0.05,  100000, 1.960}, {0.01,  100000, 2.576}, {0.001,  100000, 3.291},
        /*n == 1000000*/{0.1, 1000000, 1.645}, {0.05, 1000000, 1.960}, {0.01, 1000000, 2.576}, {0.001, 1000000, 3.291}
    };

    static const Distribution_t Khi_low_tab[] = {
        /* n == 10 */     {0.1,      10, 3.325}, {0.05,     10, 2.700}, {0.01,     10, 1.735}, {0.001,      10, 0.959},
        /* n == 100 */    {0.1,     100, 77.05}, {0.05,     100, 73.36}, {0.01,     100, 66.51}, {0.001,      100, 59.20},
        /* n == 1000 */   {0.1,    1000, 925.2}, {0.05,     1000, 913.3}, {0.01,     1000, 890.5}, {0.001,      1000, 864.9},
        /* n > 1000 Fisher approximation*/
        /* n == 10000*/  {0.1,    10000, 9768},  {0.05,  10000, 9724}, {0.01, 10000, 9639},    {0.001, 10000, 9540},
        /* n == 100000*/ {0.1,   100000, 99265},{0.05, 100000, 99125},  {0.01, 100000, 98858},  {0.001, 100000, 98544},
        /*n == 1000000*/ {0.1, 1000000, 997673}, {0.05, 1000000, 997229}, {0.01, 1000000, 996387}, {0.001, 1000000, 995394}
    };

    static const Distribution_t Khi_high_tab[] = {
        /* n == 10 */   {0.1, 10, 16.92},   {0.05, 10, 19.02},   {0.01, 10, 23.59},   {0.001, 10, 29.67},
        /* n == 100 */  {0.1, 100, 123.2},  {0.05, 100, 128.4},  {0.01, 100, 139.0},  {0.001, 100, 151.0},
        /* n == 1000 */ {0.1, 1000, 1076},  {0.05, 1000, 1089},  {0.01, 1000, 1113},  {0.001, 1000, 1142},
        /* n > 1000 */
        /* n == 10000*/{0.1, 10000, 10234},   {0.05, 10000, 10280},   {0.01, 10000, 10368},   {0.001, 10000, 10472},
        /* n == 100000*/{0.1, 100000, 100737}, {0.05, 100000, 100880}, {0.01, 100000, 101153}, {0.001, 100000, 101476},
        /* n == 1000000*/{0.1, 1000000, 1002330}, {0.05, 1000000, 1002778}, {0.01, 1000000, 1003635}, {0.001, 1000000, 1004655}
    };

    const Distribution_t* Chosen_distribution = NULL;
    size_t s = 0;
    switch(d){
        case Student: 
            Chosen_distribution = Student_tab; 
            s = *(&Student_tab + 1) - Student_tab;
            break;
        case Khi_low: 
            Chosen_distribution = Khi_low_tab; 
            s = *(&Khi_low_tab + 1) - Khi_low_tab;
            break;
        case Khi_high: 
            Chosen_distribution = Khi_high_tab; 
            s = *(&Khi_high_tab + 1) - Khi_high_tab;
            break;
        default:
        return 0.0;
    }

    //XLOOKUP Pro Max
    for(size_t i = 0; i < s; i++){
        if(Chosen_distribution[i].n == n && fabs(Chosen_distribution[i].alpha - alpha) < 1e-6) return Chosen_distribution[i].value;
    }
}

void Confidence_Interval_Var1(stat_sample_t* sample, double gamma, double* lower_limit, double *upper_limit){
    double mean = Sample_Mean(sample);
    double dev = sqrt(Sample_Variance(sample, 1));                          //bool Param = TRUE - we use unbiased variance

    double q1 = Distribution_Value(sample->size, gamma, Student);

    double q2 = (q1/(sqrt(sample->size - 1)))*dev;
    *lower_limit = mean - q2;
    *upper_limit = mean + q2;
}

void Confidence_Interval_Var2(stat_sample_t* sample, double gamma, double* lower_limit, double* upper_limit){
    double var = Sample_Variance(sample, true);

    double z1 = Distribution_Value(sample->size, gamma, Khi_low);
    double z2 = Distribution_Value(sample->size, gamma, Khi_high);

    *lower_limit = (((double)sample->size * var) / z2);
    *upper_limit = (((double)sample->size * var) / z1);
}

void Confidence_Interval_Var3(stat_sample_t* sample, double gamma, double* lower_limit, double* upper_limit){
    
    double z1 = Distribution_Value(sample->size, gamma, Normal);

    double mean, variance, st_dev;
    Sample_Analyze(sample, &mean, &variance, true, &st_dev);

    *lower_limit = mean - ((z1 * st_dev) / sqrt(sample->size));
    *upper_limit = mean + ((z1 * st_dev) / sqrt(sample->size));
}

convergence_result_t Convergence_Test_V1(uint64_t sample_size, Weibull_Type_t variation){

    uint64_t max_i_mc = 50000000;
    double z_gamma = 2.575;
    double epsilon = 0.01;
    uint64_t stabilization_n = 10000;
    uint64_t counter = 0;
    double s_q = 0;
    double s_q_sq = 0;
    double q_est, var;
    double n_star = UINT64_MAX;

    while(true){
        counter += 1;
        double current_eta = Set_Number_To_Weibull(get_rand_01(), variation);
        double sum_ksi = 0;

        for(int i = 0; i < sample_size; i++){
            double ksi_i = Set_Number_To_Weibull(get_rand_01(), Weibull_Ksi);
            sum_ksi += ksi_i;
        }

        double q_i = (sum_ksi < current_eta) ? (double)1 : (double)0;
        s_q += q_i;
        s_q_sq += q_i*q_i;

        if(counter >= 10000){
            q_est = s_q / (double)counter;
            var = ((double)1 / (counter - 1)) * (s_q_sq - counter * q_est * q_est);
            if(q_est > 0) n_star = (z_gamma * z_gamma * var) / (epsilon * epsilon * q_est* q_est);
        }

        if (counter > n_star) break;
        if (counter >= max_i_mc){
            printf("\nMonte Carlo method reached limit of iterations, params n = %d, gamma = %.2lf, execution stoped, results saved", sample_size, epsilon);
        }
    }
    convergence_result_t res;
    res.final_N = counter;
    res.variance = var;
    res.Q = q_est;
    double delta = z_gamma * sqrt(var) / sqrt((double)counter);
    res.cnf.lower_limit = q_est - delta;
    res.cnf.upper_limit = q_est + delta;

    return  res;
}

//NULL