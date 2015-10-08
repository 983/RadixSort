#include <time.h>
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <pthread.h>

#define N_THREADS 4

#ifdef __linux__
double sec(){
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9*t.tv_nsec;
}
#elif _WIN32
#include <windows.h>
double sec(){
    LARGE_INTEGER frequency, t;
    QueryPerformanceCounter(&t);
    QueryPerformanceFrequency(&frequency);
    return t.QuadPart / (double)frequency.QuadPart;
}
#else
sorry
#endif

// fast random number generator for 32 bit numbers except 0, but who needs 0 anyway
uint32_t rand32(){
    // seed
    static uint32_t x = 0x12345678;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return x;
}

// function arguments have to be passsed as a pointer for pthreads
struct radix_args {
    pthread_t thread;
    uint32_t *a, *b;
    size_t shift, n, count[256], index[256];
};

// put values into buckets based on indices
void* put_into_buckets(void *ptr){
    struct radix_args *arg = (struct radix_args*)ptr;
    uint32_t *a = arg->a, *b = arg->b;
    size_t i, shift = arg->shift, n = arg->n, *index = arg->index;

    for (i = 0; i < n; i++){
        size_t bucket = (a[i] >> shift) & 0xff;
        b[index[bucket]++] = a[i];
    }

    return NULL;
}

void* sort_range(void *ptr){
    struct radix_args *arg = (struct radix_args*)ptr;
    uint32_t *a = arg->a;
    size_t next_index, shift = arg->shift, i, n = arg->n, *count = arg->count, *index = arg->index;

    for (i = 0; i < 256; i++) count[i] = 0;

    // count numbers in buckets
    for (i = 0; i < n; i++){
        size_t bucket = (a[i] >> shift) & 0xff;
        count[bucket]++;
    }

    // accumulate indices
    next_index = 0;
    for (i = 0; i < 256; i++){
        index[i] = next_index;
        next_index += count[i];
    }

    put_into_buckets(arg);

    return NULL;
}

void radix_sort(uint32_t *a, uint32_t *b, size_t n){
    size_t i, j, shift, next_index;
    struct radix_args args[N_THREADS];
    // numbers of elements to be sorted must be multiple of number of threads
    assert(n % N_THREADS == 0);

    for (shift = 0; shift < 32; shift += 8){
        // divide elements to be sorted into N_THREAD many ranges and sort those
        for (j = 0; j < N_THREADS; j++){
            struct radix_args *arg = &args[j];
            arg->n = n / N_THREADS;
            arg->a = a + j*arg->n;
            arg->b = b + j*arg->n;
            arg->shift = shift;
            pthread_create(&arg->thread, NULL, sort_range, arg);
        }

        // wait until all threads are done
        for (j = 0; j < N_THREADS; j++){
            pthread_join(args[j].thread, NULL);
        }

        // calculate indices for buckets after merge
        next_index = 0;
        for (i = 0; i < 256; i++){
            for (j = 0; j < N_THREADS; j++){
                args[j].index[i] = next_index;
                next_index += args[j].count[i];
            }
        }

        // merge buckets of threads
        for (j = 0; j < N_THREADS; j++){
            struct radix_args *arg = &args[j];
            arg->n = n / N_THREADS;
            arg->a = b + j*arg->n;
            arg->b = a;
            arg->shift = shift;
            pthread_create(&arg->thread, NULL, put_into_buckets, arg);
        }

        // wait until all threads are done
        for (j = 0; j < N_THREADS; j++){
            pthread_join(args[j].thread, NULL);
        }
    }
}

// replace uint32_t with int to get valid Java code
void merge_sort(uint32_t vals[], uint32_t temp[], int offset, int n){
    // sorting zero or one element is easy
    if (n <= 1) return;
    // subdivide vals into left and right half
    int k, middle = n / 2;
    int i         = offset;
    int i_end     = offset + middle;
    int j         = offset + middle;
    int j_end     = offset + n;
    // sort left half
    merge_sort(vals, temp, offset, middle);
    // sort right half
    merge_sort(vals, temp, offset + middle, n - middle);
    // merge halfs into temp array until one is fully merged
    for (k = 0; i < i_end && j < j_end; k++){
        temp[k] = vals[i] <= vals[j] ? vals[i++] : vals[j++];
    }
    // merge left half into temp array if it was not merged yet
    for (; i < i_end; i++, k++) temp[k] = vals[i];
    // copy merged values from temp array back to vals
    for (i = offset; i < j; i++) vals[i] = temp[i - offset];
}

int main(){
    double t[3];
    size_t runs, i, n = 16*1024*1024;
    uint32_t *a    = (uint32_t*)malloc(sizeof(*a   )*n);
    uint32_t *b    = (uint32_t*)malloc(sizeof(*b   )*n);
    uint32_t *temp = (uint32_t*)malloc(sizeof(*temp)*n);

    for (runs = 0; runs < 5; runs++){
        // do this a few times to make sure running time is consistent

        for (i = 0; i < n; i++){
            uint32_t x = rand32();
            a[i] = x;
            b[i] = x;
        }

        t[0] = sec();

        merge_sort(a, temp, 0, n);

        t[1] = sec();

        radix_sort(b, temp, n);

        t[2] = sec();

        printf("merge_sort: %f seconds\n", t[1] - t[0]);
        printf("radix_sort: %f seconds\n", t[2] - t[1]);
        printf("\n");

        // make sure sorting algorithms are correct
        for (i = 0; i < n; i++){
            assert(a[i] == b[i]);
        }
    }

    free(a);
    free(b);
    free(temp);
    return 0;
}

