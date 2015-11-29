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

/* fast random number generator for 32 bit numbers except 0, but who needs 0 anyway */
uint32_t rand32(){
    /* seed */
    static uint32_t x = 0x12345678;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return x;
}

void radix_sort_pass(uint32_t *src, uint32_t *dst, size_t n, size_t shift){
    size_t i, next_index = 0, index[256] = {0};
    for (i = 0; i < n; i++) index[(src[i] >> shift) & 0xff]++;
    for (i = 0; i < 256; i++){
        size_t count = index[i];
        index[i] = next_index;
        next_index += count;
    }
    for (i = 0; i < n; i++) dst[index[(src[i] >> shift) & 0xff]++] = src[i];
}

void radix_sort(uint32_t *a, uint32_t *temp, size_t n){
    radix_sort_pass(a, temp, n, 0*8);
    radix_sort_pass(temp, a, n, 1*8);
    radix_sort_pass(a, temp, n, 2*8);
    radix_sort_pass(temp, a, n, 3*8);
}

/* function arguments have to be passsed as a pointer for pthreads */
struct radix_args {
    pthread_t thread;
    uint32_t *a, *b;
    size_t shift, n, count[256], index[256];
};

/* put values into buckets based on indices */
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

    /* count numbers in buckets */
    for (i = 0; i < n; i++){
        size_t bucket = (a[i] >> shift) & 0xff;
        count[bucket]++;
    }

    /* accumulate indices */
    next_index = 0;
    for (i = 0; i < 256; i++){
        index[i] = next_index;
        next_index += count[i];
    }

    put_into_buckets(arg);

    return NULL;
}

void parallel_radix_sort(uint32_t *a, uint32_t *b, size_t n){
    size_t i, j, shift, next_index;
    struct radix_args args[N_THREADS];
    /* numbers of elements to be sorted must be multiple of number of threads */
    assert(n % N_THREADS == 0);

    for (shift = 0; shift < 32; shift += 8){
        /* divide elements to be sorted into N_THREAD many ranges and sort those */
        for (j = 0; j < N_THREADS; j++){
            struct radix_args *arg = &args[j];
            arg->n = n / N_THREADS;
            arg->a = a + j*arg->n;
            arg->b = b + j*arg->n;
            arg->shift = shift;
            pthread_create(&arg->thread, NULL, sort_range, arg);
        }

        /* wait until all threads are done */
        for (j = 0; j < N_THREADS; j++){
            pthread_join(args[j].thread, NULL);
        }

        /* calculate indices for buckets after merge */
        next_index = 0;
        for (i = 0; i < 256; i++){
            for (j = 0; j < N_THREADS; j++){
                args[j].index[i] = next_index;
                next_index += args[j].count[i];
            }
        }

        /* merge buckets of threads */
        for (j = 0; j < N_THREADS; j++){
            struct radix_args *arg = &args[j];
            arg->n = n / N_THREADS;
            arg->a = b + j*arg->n;
            arg->b = a;
            arg->shift = shift;
            pthread_create(&arg->thread, NULL, put_into_buckets, arg);
        }

        /* wait until all threads are done */
        for (j = 0; j < N_THREADS; j++){
            pthread_join(args[j].thread, NULL);
        }
    }
}

void merge_sort(uint32_t vals[], uint32_t temp[], int offset, int n){
    int k, middle, i, i_end, j, j_end;
    /* sorting zero or one element is easy */
    if (n <= 1) return;
    /* subdivide vals into left and right half */
    middle    = n / 2;
    i         = offset;
    i_end     = offset + middle;
    j         = offset + middle;
    j_end     = offset + n;
    /* sort left half */
    merge_sort(vals, temp, offset, middle);
    /* sort right half */
    merge_sort(vals, temp, offset + middle, n - middle);
    /* merge halfs into temp array until one is fully merged */
    for (k = 0; i < i_end && j < j_end; k++){
        temp[k] = vals[i] <= vals[j] ? vals[i++] : vals[j++];
    }
    /* merge left half into temp array if it was not merged yet */
    for (; i < i_end; i++, k++) temp[k] = vals[i];
    /* copy merged values from temp array back to vals */
    for (i = offset; i < j; i++) vals[i] = temp[i - offset];
}

int compare_uint32_t(const void *a, const void *b){
    uint32_t x = *(uint32_t*)a;
    uint32_t y = *(uint32_t*)b;
    return x == y ? 0 : (x < y ? -1 : +1);
}

uint32_t *a, *b, *c, *temp;

size_t n;

void init_random(){
    size_t i;
    for (i = 0; i < n; i++){
        uint32_t x = rand32();
        a[i] = x;
        b[i] = x;
        c[i] = x;
    }
}

void init_sorted(){
    size_t i;
    for (i = 0; i < n; i++){
        a[i] = i;
        b[i] = i;
        c[i] = i;
    }
}

void init_reversed(){
    size_t i;
    for (i = 0; i < n; i++){
        size_t j = n - 1 - i;
        a[i] = j;
        b[i] = j;
        c[i] = j;
    }
}

void init_90_percent_sorted(){
    size_t i;
    init_sorted();
    for (i = 0; i < n/10; i++){
        uint32_t x = rand32();
        size_t j = rand32() % n;
        a[j] = x;
        b[j] = x;
        c[j] = x;
    }
}

#define N_ALGORITHMS 4

void benchmark(void (*init_func)()){
    printf("n, qsort, merge_sort, radix_sort, parallel_radix_sort\n");
    for (n = N_THREADS; n <= 4*1024*1024; n *= 2){
        double dt, t[N_ALGORITHMS + 1];
        size_t i, runs;
        double fastest[N_ALGORITHMS];
        for (i = 0; i < N_ALGORITHMS; i++) fastest[i] = 1e20;

        a    = (uint32_t*)malloc(sizeof(*a   )*n);
        b    = (uint32_t*)malloc(sizeof(*b   )*n);
        c    = (uint32_t*)malloc(sizeof(*c   )*n);
        temp = (uint32_t*)malloc(sizeof(*temp)*n);

         /* do this a few times to make sure running time is consistent */
        for (runs = 0; runs < 10; runs++){

            init_func();

            t[0] = sec();

            qsort(a, n, sizeof(uint32_t), compare_uint32_t);

            t[1] = sec();

            merge_sort(a, temp, 0, n);

            t[2] = sec();

            radix_sort(b, temp, n);

            t[3] = sec();

            parallel_radix_sort(c, temp, n);

            t[4] = sec();

            for (i = 0; i < N_ALGORITHMS; i++){
                dt = t[i + 1] - t[i];
                if (dt < fastest[i]) fastest[i] = dt;
            }

            /* make sure sorting algorithms are correct */
            for (i = 0; i < n; i++){
                assert(a[i] == b[i]);
                assert(a[i] == c[i]);
            }
        }

        printf("%9i, ", (int)n);
        for (i = 0; i < N_ALGORITHMS; i++) printf("%f, ", fastest[i]);
        printf("\n");

        free(a);
        free(b);
        free(c);
        free(temp);
    }
}

int main(){
    printf("sorted\n\n");
    benchmark(init_sorted);
    printf("\nreversed\n\n");
    benchmark(init_reversed);
    printf("\n90 percent sorted\n\n");
    benchmark(init_90_percent_sorted);
    printf("\nrandom\n\n");
    benchmark(init_random);
    return 0;
}

