#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include <x86intrin.h>

#define PI 3.14159265358979323846

static void
bmark(void (func)(float *amp, int K, float *buff, int N), char *name, int iters)
{
    printf("Function: %s\n", name);
    int K = 1200;
    int N = 2048;

    float amp[K];
    float buff[N];

    for (int i = 0; i < K; i++) {
        amp[i] = 1.0f/K;
    }
    amp[0] = 0.0f;

    // memset(buff, 0, sizeof(float)*N);

    clock_t t0 = clock();
    for (int i = 0; i < iters; i++) {
        func(amp, K, buff, N);
    }
    clock_t t1 = clock();

    double total_time = ((double)(t1-t0))/(CLOCKS_PER_SEC);
    double avg_time = total_time/iters;
    printf("K=%d, N=%d, total: %0.15lf s avg: %0.15lf s\n", K, N, total_time, avg_time);

    for (int i = 0; i < N; i++) {
        fprintf(stderr, "%f ", buff[i]);
    }
    fprintf(stderr, "\n");

    printf("\n");
}

static void
direct(float *amp, int K, float *buff, int N)
{
    float t = 0.0f;
    float p = 1.0f/48000.0f;
    for (int n = 0; n < N; n++) {
        buff[n] = 0.0f;
//        printf("n: %d\t", n);
        for (int k = 0; k < K; k++) {
            float v = sinf(k*(t+n*p));
            buff[n] += amp[k]*v;
//            printf("\t(k: %d, %f) ", k, v);
        }
///        printf("\n");
    }
}

/*
    Sum a sine table, iterating through it at different speeds for different
    freqs.

    O(LK), where L is length of sine table?
*/
float *G_SINETABLE;
int G_SINETABLE_LEN;
static void
wave_table(float *amp, int K, float *buff, int N)
{
    for (int k = 0; k < K; k++) {
        int i = 0;
        for (int n = 0; n < N; n++) {
            /* no interpolation... yet */
            buff[n] += amp[k]*G_SINETABLE[i];
            i = (i + k*G_SINETABLE_LEN/48000) % G_SINETABLE_LEN;
        }
    }
}

/*
    Inverse FFT

    O(NlogN)
*/
void
ifft(float *amp, int K, float *buff, int N)
{
}

/*
    Compute using additive and angle multiple trig identities.
*/
static void
m2(float *amp, int K, float *buff, int N)
{
    float t = 0.0f;
    float p = 1.0f/48000.0f;

    float sin_p = sinf(p);
    float cos_p = cosf(p);
    float sin_t = 0.0f;
    float cos_t = 1.0f;

    float sin_1, sin_2, sin_tmp, n_cos_t, n_sin_t;

    for (int n = 0; n < N; n++) {
        t = p*n;
        sin_2 = 0.0f;
        sin_1 = sin_t;

        buff[n] = 0.0f;
        for (int k = 0; k < K; k++) {
            sin_tmp = 2.0f*sin_1*cos_t - sin_2;
            buff[n] += sin_tmp;
            sin_2 = sin_1;
            sin_1 = sin_tmp;
        }

        n_sin_t = sin_t*cos_p + sin_p*cos_t;
        n_cos_t = cos_t*cos_p - sin_t*sin_p;
        cos_t = n_cos_t;
        sin_t = n_sin_t;
    }
}

static void
m2_par2(float *amp, int K, float *buff, int N)
{

    float t = 0.0f; // t is the time of the start of the buffer.
    float p = 1.0f/48000.0f;

    // Compute harmonics tables w/ SSE
    __m128 skt[K/4];
    __m128 ckt[K/4];

    float s_skt_0[4];
    float s_ckt_0[4];
    for (int i = 0; i < 4; i++) {
        s_skt_0[i] = sinf(i*t);
        s_ckt_0[i] = cosf(i*t);
    }    

    float s_s4t = sinf(4.0f*t);
    float s_c4t = cosf(4.0f*t);

    __m128 s4t = _mm_load1_ps(&s_s4t);
    __m128 c4t = _mm_load1_ps(&s_c4t);

    skt[0] = _mm_loadu_ps(s_skt_0);
    ckt[0] = _mm_loadu_ps(s_ckt_0);

    __m128 skt_s4t, ckt_c4t, skt_c4t, ckt_s4t;
    for (int k4 = 1; k4 < K/4; k4++) {
        skt_s4t = _mm_mul_ps(skt[k4-1], s4t);
        skt_c4t = _mm_mul_ps(skt[k4-1], c4t);
        ckt_s4t = _mm_mul_ps(ckt[k4-1], s4t);
        ckt_c4t = _mm_mul_ps(ckt[k4-1], c4t);

        skt[k4] = _mm_add_ps(skt_c4t, ckt_s4t);
        ckt[k4] = _mm_sub_ps(ckt_c4t, skt_s4t);
    }


    // Each output sample is the sum of harmonics at the time, scaled by amp
    for (int n = 0; n < N; n++) {
        buff[n] = 0.0f;
        for (int k4 = 0; k4 < K/4; k4++) {
            float tmp_sink[4];
            float tmp_cosk[4];
            _mm_store_ps(tmp_sink, skt[k4]);
            _mm_store_ps(tmp_cosk, ckt[k4]);

            for (int i = 0; i < 4; i++) {
                int k = i*4;
                // buff[n] += amp[k]*(tmp_sink[i]*cosf(k*n*p) + tmp_cosk[i]*sinf(k*n*p));
                buff[n] += amp[k]*(sinf(k*t)*cosf(k*n*p) + cosf(k*t)*sinf(k*n*p));
            }
        }
    }
}

static void
m4(float *amp, int K, float* buff, int N)
{

    /* Using double instead of float... Unfortunately using float itself seemed to introduce more error. */
    float t = 0.0;
    float p = 1.0/48000.0;
    float sin_knp, cos_knp, sin_np, cos_np, sin_kt, cos_kt, sin_p, cos_p, sin_t, cos_t, sin_kp, cos_kp;
    float n_sin_kt, n_cos_kt, n_sin_knp, n_cos_knp, n_sin_kp, n_cos_kp, n_sin_np, n_cos_np;

    sin_np = 0.0; // sinf(0*p);
    cos_np = 1.0; // cosf(0*p);

    sin_p = sin(p);
    cos_p = cos(p);
    sin_t = sin(t);
    cos_t = cos(t);

    float buff_n;

    for (int n = 0; n < N; n++) {
        buff_n = 0.0f;
        sin_kt = 0.0f; // sinf(0*t);
        cos_kt = 1.0f; // cosf(0*t);
        sin_knp = 0.0f;
        cos_knp = 1.0f;
        sin_kp = 0.0f;
        cos_kp = 1.0f;

        for (int k = 0; k < K; k++) {
            buff_n += amp[k]*(sin_kt*cos_knp + cos_kt*sin_knp);
            n_sin_kt = sin_kt*cos_t + cos_kt*sin_t;
            n_cos_kt = cos_kt*cos_t - sin_kt*sin_t;
            n_sin_knp = sin_knp*cos_np + cos_knp*sin_np;
            n_cos_knp = cos_knp*cos_np - sin_knp*sin_np;
            n_sin_kp = sin_kp*cos_p + cos_kp*sin_p;
            n_cos_kp = cos_kp*cos_p - sin_kp*sin_p;

            sin_kt = n_sin_kt;
            cos_kt = n_cos_kt;
            sin_kp = n_sin_kp;
            cos_kp = n_cos_kp;
            sin_knp = n_sin_knp;
            cos_knp = n_cos_knp;
        }

        n_sin_np = sin_np*cos_p + cos_np*sin_p;
        n_cos_np = cos_np*cos_p - sin_np*sin_p;
        sin_np = n_sin_np;
        cos_np = n_cos_np;

        buff[n] = buff_n;
    }
}


static void
m4_arr(float *amp, int K, float* buff, int N)
{

    /* Using double instead of float... Unfortunately using float itself seemed to introduce more error. */
    float t = 0.0;
    float p = 1.0/48000.0;
    float sin_p, cos_p, sin_t, cos_t;
    float sin_np[N];
    float cos_np[N];

    float sin_kt[K];
    float cos_kt[K];

    float sin_knp[K];
    float cos_knp[K];

    float n_sin_knp, n_cos_knp;

    sin_np[0] = 0.0; // sinf(0*p);
    cos_np[0] = 1.0; // cosf(0*p);

    sin_p = sin(p);
    cos_p = cos(p);
    sin_t = sin(t);
    cos_t = cos(t);

    float buff_n;

    for (int n = 1; n < N; n++) {
        sin_np[n] = sin_np[n-1]*cos_p + cos_np[n-1]*sin_p;
        cos_np[n] = cos_np[n-1]*cos_p - sin_np[n-1]*sin_p;
    }

    sin_kt[0] = 0.0f; // sinf(0*t);
    cos_kt[0] = 1.0f; // cosf(0*t);

    for (int k = 1; k < K; k++) {
        sin_kt[k] = sin_kt[k-1]*cos_t + cos_kt[k-1]*sin_t;
        cos_kt[k] = cos_kt[k-1]*cos_t - sin_kt[k-1]*sin_t;
    }

    sin_knp[0] = 0.0f;
    cos_knp[0] = 1.0f;
    for (int n = 0; n < N; n++) {
        for (int k = 1; k < K; k++) {
            sin_knp[k] = sin_knp[k-1]*cos_np[n] + cos_knp[k-1]*sin_np[n];
            cos_knp[k] = cos_knp[k-1]*cos_np[n] - sin_knp[k-1]*sin_np[n];
        }

        buff[n] = 0.0f;

        for (int k = 0; k < K; k++) {
            buff[n] += amp[k]*(sin_kt[k]*cos_knp[k] + cos_kt[k]*sin_knp[k]);
        }
    }
}


static void
m4_sse(float *amp, int K, float* buff, int N)
{

    float t = 0.0;
    float p = 1.0/48000.0;
    float sin_p, cos_p, sin_t, cos_t;

    float sin_np[N];
    float cos_np[N];

    sin_p = sin(p);
    cos_p = cos(p);
    sin_t = sin(t);
    cos_t = cos(t);
    __m128 sse_cos_p = _mm_set1_ps(cos_p);
    __m128 sse_sin_p = _mm_set1_ps(sin_p);
    __m128 sse_cos_t = _mm_set1_ps(cos_t);
    __m128 sse_sin_t = _mm_set1_ps(sin_t);

    float buff_n;

    __m128 sse_sin_np[N/4];
    __m128 sse_cos_np[N/4];

    sse_sin_np[0] = _mm_setr_ps(0.0f, sin_p, sinf(2*p), sinf(3*p));
    sse_cos_np[0] = _mm_setr_ps(1.0f, cos_p, cosf(2*p), cosf(3*p));
    _mm_store_ps(sin_np, sse_sin_np[0]);
    _mm_store_ps(cos_np, sse_cos_np[0]);

    for (int n4 = 1; n4 < N/4; n4++) {
        __m128 snp_cp = _mm_mul_ps(sse_sin_np[n4-1], sse_cos_p);
        __m128 snp_sp = _mm_mul_ps(sse_sin_np[n4-1], sse_sin_p);
        __m128 cnp_cp = _mm_mul_ps(sse_cos_np[n4-1], sse_cos_p);
        __m128 cnp_sp = _mm_mul_ps(sse_cos_np[n4-1], sse_sin_p);

        sse_sin_np[n4] = _mm_add_ps(snp_cp, cnp_sp);
        sse_cos_np[n4] = _mm_sub_ps(cnp_cp, snp_sp);

        _mm_store_ps(sin_np + 4*n4, sse_sin_np[n4]);
        _mm_store_ps(cos_np + 4*n4, sse_cos_np[n4]);
    }

    __m128 sse_sin_kt[K/4];
    __m128 sse_cos_kt[K/4];

    sse_sin_kt[0] = _mm_setr_ps(0.0f, sin_t, sinf(2*t), sinf(3*t));
    sse_cos_kt[0] = _mm_setr_ps(1.0f, cos_t, cosf(2*t), cosf(3*t));

    for (int k4 = 1; k4 < K/4; k4++) {
        __m128 skt_ct = _mm_mul_ps(sse_sin_kt[k4-1], sse_cos_t);
        __m128 skt_st = _mm_mul_ps(sse_sin_kt[k4-1], sse_sin_t);
        __m128 ckt_ct = _mm_mul_ps(sse_cos_kt[k4-1], sse_cos_t);
        __m128 ckt_st = _mm_mul_ps(sse_cos_kt[k4-1], sse_sin_t);

        sse_sin_kt[k4] = _mm_add_ps(skt_ct, ckt_st);
        sse_cos_kt[k4] = _mm_sub_ps(ckt_ct, skt_st);
    }

    __m128 sse_sin_knp[K/4];
    __m128 sse_cos_knp[K/4];
    __m128 sse_sin_ktnp[K/4];

    float sin_ktnp[K];

    for (int n = 0; n < N; n++) {
        sse_sin_knp[0] = _mm_setr_ps(0.0f, sinf(n*p), sinf(2.0f*n*p), sinf(3.0f*n*p));
        sse_cos_knp[0] = _mm_setr_ps(1.0f, cosf(n*p), cosf(2.0f*n*p), cosf(3.0f*n*p));

        __m128 curr_cos_np = _mm_set1_ps(cos_np[n]);
        __m128 curr_sin_np = _mm_set1_ps(sin_np[n]);

        {
            __m128 skt_cknp = _mm_mul_ps(sse_sin_kt[0], sse_cos_knp[0]);
            __m128 ckt_sknp = _mm_mul_ps(sse_cos_kt[0], sse_sin_knp[0]);

            sse_sin_ktnp[0] = _mm_add_ps(skt_cknp, ckt_sknp);

            _mm_store_ps(sin_ktnp, sse_sin_ktnp[0]);
        }

        for (int k4 = 1; k4 < K/4; k4++) {
            __m128 sknp_cnp = _mm_mul_ps(sse_sin_knp[k4-1], curr_cos_np);
            __m128 sknp_snp = _mm_mul_ps(sse_sin_knp[k4-1], curr_sin_np);
            __m128 cknp_cnp = _mm_mul_ps(sse_cos_knp[k4-1], curr_cos_np);
            __m128 cknp_snp = _mm_mul_ps(sse_cos_knp[k4-1], curr_sin_np);

            sse_sin_knp[k4] = _mm_add_ps(sknp_cnp, cknp_snp);
            sse_cos_knp[k4] = _mm_sub_ps(cknp_cnp, sknp_snp);

            __m128 skt_cknp = _mm_mul_ps(sse_sin_kt[k4], sse_cos_knp[k4]);
            __m128 ckt_sknp = _mm_mul_ps(sse_cos_kt[k4], sse_sin_knp[k4]);

            sse_sin_ktnp[k4] = _mm_add_ps(skt_cknp, ckt_sknp);

            _mm_store_ps(sin_ktnp + 4*k4, sse_sin_ktnp[k4]);
        }

        buff[n] = 0.0f;
        for (int k = 0; k < K; k++) {
            buff[n] += amp[k]*sin_ktnp[k];
        }
    }
}

static void
m4_sse_opt(float *amp, int K, float* buff, int N)
{

    // Vars
    float t = 0.0;
    float p = 1.0/48000.0;
    float sin_p = sin(p);
    float cos_p = cos(p);
    float sin_t = sin(t);
    float cos_t = cos(t);
    __m128 sse_cos_p = _mm_set1_ps(cos_p);
    __m128 sse_sin_p = _mm_set1_ps(sin_p);
    __m128 sse_cos_t = _mm_set1_ps(cos_t);
    __m128 sse_sin_t = _mm_set1_ps(sin_t);

    // Arrays
    float sin_np[N];
    float cos_np[N];
    __m128 sse_sin_kt[K/4];
    __m128 sse_cos_kt[K/4];
 //   __m128 sse_sin_np[N/4];
//    __m128 sse_cos_np[N/4];
    __m128 sse_amp[K/4];
    
    for (int k4 = 0; k4 < K/4; k4++) {
        sse_amp[k4] = _mm_load_ps(amp+4*k4);
    }

/*
    sse_sin_np[0] = _mm_setr_ps(0.0f, sin_p, sinf(2*p), sinf(3*p));
    sse_cos_np[0] = _mm_setr_ps(1.0f, cos_p, cosf(2*p), cosf(3*p));
    _mm_store_ps(sin_np, sse_sin_np[0]);
    _mm_store_ps(cos_np, sse_cos_np[0]);

    for (int n4 = 1; n4 < N/4; n4++) {
        __m128 snp_cp = _mm_mul_ps(sse_sin_np[n4-1], sse_cos_p);
        __m128 snp_sp = _mm_mul_ps(sse_sin_np[n4-1], sse_sin_p);
        __m128 cnp_cp = _mm_mul_ps(sse_cos_np[n4-1], sse_cos_p);
        __m128 cnp_sp = _mm_mul_ps(sse_cos_np[n4-1], sse_sin_p);

        sse_sin_np[n4] = _mm_add_ps(snp_cp, cnp_sp);
        sse_cos_np[n4] = _mm_sub_ps(cnp_cp, snp_sp);

        _mm_store_ps(sin_np + 4*n4, sse_sin_np[n4]);
        _mm_store_ps(cos_np + 4*n4, sse_cos_np[n4]);
    }
*/
    for (int n = 0; n < N; n++) {
        sin_np[n] = sinf(n*p);
        cos_np[n] = cosf(n*p);
    }

    sse_sin_kt[0] = _mm_setr_ps(0.0f, sin_t, sinf(2*t), sinf(3*t));
    sse_cos_kt[0] = _mm_setr_ps(1.0f, cos_t, cosf(2*t), cosf(3*t));

    for (int k4 = 1; k4 < K/4; k4++) {
        __m128 skt_ct = _mm_mul_ps(sse_sin_kt[k4-1], sse_cos_t);
        __m128 skt_st = _mm_mul_ps(sse_sin_kt[k4-1], sse_sin_t);
        __m128 ckt_ct = _mm_mul_ps(sse_cos_kt[k4-1], sse_cos_t);
        __m128 ckt_st = _mm_mul_ps(sse_cos_kt[k4-1], sse_sin_t);

        // sse_sin_kt[k4] = _mm_add_ps(skt_ct, ckt_st);
        // sse_cos_kt[k4] = _mm_sub_ps(ckt_ct, skt_st);

        int k = 4*k4;
        sse_sin_kt[k4] = _mm_setr_ps(sinf(k*t), sinf((k+1)*t), sinf((k+2)*t), sinf((k+3)*t));
        sse_cos_kt[k4] = _mm_setr_ps(cosf(k*t), cosf((k+1)*t), cosf((k+2)*t), cosf((k+3)*t));
        // TODO: Change and test this with changing t
        // sse_sin_kt[k4] = _mm_set1_ps(0.0f);
        // sse_cos_kt[k4] = _mm_set1_ps(1.0f);
    }

    __m128 sse_sin_knp[K/4];
    __m128 sse_cos_knp[K/4];
    __m128 sse_sin_ktnp;

    for (int n = 0; n < N; n++) {
        sse_sin_knp[0] = _mm_setr_ps(0.0f, sinf(n*p), sinf(2.0f*n*p), sinf(3.0f*n*p));
        sse_cos_knp[0] = _mm_setr_ps(1.0f, cosf(n*p), cosf(2.0f*n*p), cosf(3.0f*n*p));

        __m128 curr_cos_np = _mm_set1_ps(cosf(n*p));
        __m128 curr_sin_np = _mm_set1_ps(sinf(n*p));
        __m128 skt_cknp = _mm_mul_ps(sse_sin_kt[0], sse_cos_knp[0]);
        __m128 ckt_sknp = _mm_mul_ps(sse_cos_kt[0], sse_sin_knp[0]);

        sse_sin_ktnp = _mm_add_ps(skt_cknp, ckt_sknp);

        buff[n] = 0.0f;
        __m128 mix_out = _mm_mul_ps(sse_amp[0], sse_sin_ktnp);
        for (int k4 = 1; k4 < K/4; k4++) {
            int k = 4*k4;

            //__m128 sknp_snp = _mm_mul_ps(sse_sin_knp[k4-1], curr_sin_np);
            //__m128 cknp_cnp = _mm_mul_ps(sse_cos_knp[k4-1], curr_cos_np);
            //__m128 cknp_snp = _mm_mul_ps(sse_cos_knp[k4-1], curr_sin_np);
            //__m128 sknp_cnp = _mm_mul_ps(sse_sin_knp[k4-1], curr_cos_np);

            float a[4];
            _mm_storeu_ps(a, sse_sin_knp[k4-1]);

            __m128 sknp_cnp = _mm_setr_ps(sinf(k*n*p), sinf((k+1)*n*p), sinf((k+2)*n*p), sinf((k+3)*n*p));
            sknp_cnp = _mm_mul_ps(sknp_cnp, curr_cos_np);

            __m128 cknp_cnp = _mm_setr_ps(cosf(k*n*p), cosf((k+1)*n*p), cosf((k+2)*n*p), cosf((k+3)*n*p));
            cknp_cnp = _mm_mul_ps(cknp_cnp, curr_cos_np);

            __m128 sknp_snp = _mm_setr_ps(sinf(k*n*p), sinf((k+1)*n*p), sinf((k+2)*n*p), sinf((k+3)*n*p));
            sknp_snp = _mm_mul_ps(sknp_cnp, curr_sin_np);

            __m128 cknp_snp = _mm_setr_ps(cosf(k*n*p), cosf((k+1)*n*p), cosf((k+2)*n*p), cosf((k+3)*n*p));
            cknp_snp = _mm_mul_ps(cknp_snp, curr_sin_np);

            sse_sin_knp[k4] = _mm_add_ps(sknp_cnp, cknp_snp);
            sse_cos_knp[k4] = _mm_sub_ps(cknp_cnp, sknp_snp);

            // sse_sin_knp[k4] = _mm_setr_ps(sinf(k*n*p), sinf((k+1)*n*p), sinf((k+2)*n*p), sinf((k+3)*n*p));

            skt_cknp = _mm_mul_ps(sse_sin_kt[k4], sse_cos_knp[k4]);
            ckt_sknp = _mm_mul_ps(sse_cos_kt[k4], sse_sin_knp[k4]);

            sse_sin_ktnp = _mm_add_ps(skt_cknp, ckt_sknp);
            __m128 a_sin_ktnp = _mm_mul_ps(sse_amp[k4], sse_sin_ktnp);

            mix_out = _mm_add_ps(mix_out, a_sin_ktnp);
        }

        // horiztonal sum of mix_out to buff[n]
        // SSE3!
        mix_out = _mm_hadd_ps(mix_out, mix_out);
        mix_out = _mm_hadd_ps(mix_out, mix_out);
        _mm_store_ss(buff+n, mix_out);
    }
}

static void
m4_sse_opt2(float *amp, int K, float* buff, int N)
{
    // Vars
    float t = 0.0;
    float p = 1.0/48000.0;
    float sin_p = sin(p);
    float cos_p = cos(p);
    float sin_t = sin(t);
    float cos_t = cos(t);
    __m128 sse_amp;
    __m128 sse_cos_p = _mm_set1_ps(cos_p);
    __m128 sse_sin_p = _mm_set1_ps(sin_p);
    __m128 sse_cos_t = _mm_set1_ps(cos_t);
    __m128 sse_sin_t = _mm_set1_ps(sin_t);
    __m128 sse_sin_ktnp;
    __m128 snp_cp, snp_sp, cnp_cp, cnp_sp;
    __m128 skt_ct, skt_st, ckt_ct, ckt_st;
    __m128 mix_out;
    __m128 sse_sin_knp_1, sse_cos_knp_1, curr_cos_np, curr_sin_np;
    __m128 sknp_cnp, sknp_snp, cknp_cnp, cknp_snp, sse_sin_knp, sse_cos_knp, skt_cknp, ckt_sknp, a_sin_ktnp;

    // Arrays
    float sin_np[N];
    float cos_np[N];
    __m128 sse_sin_kt[K/4];
    __m128 sse_cos_kt[K/4];
    // __m128 sse_amp[K/4];
    __m128 sse_sin_np[N/4];
    __m128 sse_cos_np[N/4];
/*    
    for (int k4 = 0; k4 < K/4; k4++) {
        sse_amp[k4] = _mm_load_ps(amp+4*k4);
    }
*/

    // Might need setr instead of set
    sse_sin_np[0] = _mm_setr_ps(0.0f, sin_p, sinf(2*p), sinf(3*p)); // TODO: remove sinf
    sse_cos_np[0] = _mm_setr_ps(1.0f, cos_p, cosf(2*p), cosf(3*p)); // TODO: remove cosf
    _mm_store_ps(sin_np, sse_sin_np[0]);
    _mm_store_ps(cos_np, sse_cos_np[0]);

    for (int n4 = 1; n4 < N/4; n4++) {
        snp_cp = _mm_mul_ps(sse_sin_np[n4-1], sse_cos_p);
        snp_sp = _mm_mul_ps(sse_sin_np[n4-1], sse_sin_p);
        cnp_cp = _mm_mul_ps(sse_cos_np[n4-1], sse_cos_p);
        cnp_sp = _mm_mul_ps(sse_cos_np[n4-1], sse_sin_p);

        sse_sin_np[n4] = _mm_add_ps(snp_cp, cnp_sp);
        sse_cos_np[n4] = _mm_sub_ps(cnp_cp, snp_sp);

        _mm_store_ps(sin_np + 4*n4, sse_sin_np[n4]);
        _mm_store_ps(cos_np + 4*n4, sse_cos_np[n4]);
    }

    // Might need setr instead of set
    sse_sin_kt[0] = _mm_setr_ps(0.0f, sin_t, sinf(2*t), sinf(3*t)); // TODO: remove sinf
    sse_cos_kt[0] = _mm_setr_ps(1.0f, cos_t, cosf(2*t), cosf(3*t)); // TODO: remove cosf

    for (int k4 = 1; k4 < K/4; k4++) {
        skt_ct = _mm_mul_ps(sse_sin_kt[k4-1], sse_cos_t);
        skt_st = _mm_mul_ps(sse_sin_kt[k4-1], sse_sin_t);
        ckt_ct = _mm_mul_ps(sse_cos_kt[k4-1], sse_cos_t);
        ckt_st = _mm_mul_ps(sse_cos_kt[k4-1], sse_sin_t);

        sse_sin_kt[k4] = _mm_add_ps(skt_ct, ckt_st);
        sse_cos_kt[k4] = _mm_sub_ps(ckt_ct, skt_st);
    }

    for (int n = 0; n < N; n++) {
        //Might need setr instead of set
        sse_sin_knp_1 = _mm_setr_ps(0.0f, sinf(n*p), sinf(2.0f*n*p), sinf(3.0f*n*p)); //TODO: compute these values w/ out calling sinf
        sse_cos_knp_1 = _mm_setr_ps(1.0f, cosf(n*p), cosf(2.0f*n*p), cosf(3.0f*n*p)); //TODO: compute these values w/out cosf

        curr_cos_np = _mm_set1_ps(cos_np[n]);
        curr_sin_np = _mm_set1_ps(sin_np[n]);

        skt_cknp = _mm_mul_ps(sse_sin_kt[0], sse_cos_knp_1);
        ckt_sknp = _mm_mul_ps(sse_cos_kt[0], sse_sin_knp_1);

        sse_sin_ktnp = _mm_add_ps(skt_cknp, ckt_sknp);

        buff[n] = 0.0f;

        sse_amp = _mm_load_ps(amp);
        mix_out = _mm_mul_ps(sse_amp, sse_sin_ktnp);

        for (int k4 = 1; k4 < K/4; k4++) {
            sknp_cnp = _mm_mul_ps(sse_sin_knp_1, curr_cos_np);
            sknp_snp = _mm_mul_ps(sse_sin_knp_1, curr_sin_np);
            cknp_cnp = _mm_mul_ps(sse_cos_knp_1, curr_cos_np);
            cknp_snp = _mm_mul_ps(sse_cos_knp_1, curr_sin_np);

            sse_sin_knp = _mm_add_ps(sknp_cnp, cknp_snp);
            sse_cos_knp = _mm_sub_ps(cknp_cnp, sknp_snp);

            skt_cknp = _mm_mul_ps(sse_sin_kt[k4], sse_cos_knp);
            ckt_sknp = _mm_mul_ps(sse_cos_kt[k4], sse_sin_knp);

            sse_sin_ktnp = _mm_add_ps(skt_cknp, ckt_sknp);

            sse_amp = _mm_load_ps(amp+4*k4);
            a_sin_ktnp = _mm_mul_ps(sse_amp, sse_sin_ktnp);
            
            mix_out = _mm_add_ps(mix_out, a_sin_ktnp);
        }

        // horiztonal sum of mix_out to buff[n]
        // SSE3!
        mix_out = _mm_hadd_ps(mix_out, mix_out);
        mix_out = _mm_hadd_ps(mix_out, mix_out);
        _mm_store_ss(buff+n, mix_out);
    }
}
/*
    (k+1)*x = kx + x

    w = 2*pi*f
    w*(t+(n+1)*p) = w*(t + np + p) = w*(t + np) + w*p

    w2 = w+d

    w2(t+np) = (w+d)(t+np) = w(t+np) + d(t+np)
*/
static void
m4_r(float *amp, int K, float *buff, int N)
{
    float p = 1.0f/48000.0f;
    float t = 0.0f;
    float f = 1.0f/(2.0f*PI);
    float w = 2 * PI * f;
    
    float swp = sinf(w*p);
    float cwp = cosf(w*p);
    float st = sinf(t);
    float ct = cosf(t);
    float ckx[K];
    float skx[K];
    float cx[N];
    float sx[N];


    // TODO: change this to work when frequency changes
    // To change frequency, simply update cwp and swp in the loop
    // cwp = cosf(w2), swp = sinf(w2)
    cx[0] = ct;
    sx[0] = st;
    for (int n = 1; n < N; n++) {
        // float x = w*(t+n*p);
        cx[n] = cx[n-1]*cwp - sx[n-1]*swp;
        sx[n] = sx[n-1]*cwp + cx[n-1]*swp;
    }

    for (int n = 0; n < N; n++) {

        skx[0] = 0.0f;
        ckx[0] = 1.0f;
        buff[n] = amp[0]*skx[0];
        for (int k = 1; k < K; k++) {
            skx[k] = skx[k-1]*cx[n] + ckx[k-1]*sx[n];
            ckx[k] = ckx[k-1]*cx[n] - skx[k-1]*sx[n];

            buff[n] += amp[k]*skx[k]; 
        }
    }
}

static void
m4_rs(float *amp, int K, float *buff, int N)
{
    float p, t, f, w, swp, cwp, st, ct, ckx, ckx_1, skx, skx_1, cx, cx_1, sx, sx_1, sum;

    p = 1.0f/48000.0f;
    t = 0.0f;
    f = 1.0f/(2.0f*PI);
    w = 2 * PI * f;
    
    swp = sinf(w*p);
    cwp = cosf(w*p);
    st = sinf(t);
    ct = cosf(t);

    // TODO: change this to work when frequency changes
    // To change frequency, simply update cwp and swp in the loop
    // cwp = cosf(w2), swp = sinf(w2)
    // TODO: suppot changing amp[k] on a sample-by-sample basis
    // Common case:
    // Changing f often
    // Changing amp[k] not so much.

    cx = ct;
    sx = st;
    for (int n = 0; n < N; n++) {
        skx_1 = 0.0f;
        ckx_1 = 1.0f;
        sum = amp[0]*skx_1;

        for (int k = 1; k < K; k++) {
            skx = skx_1*cx + ckx_1*sx;
            ckx = ckx_1*cx - skx_1*sx;
            skx_1 = skx;
            ckx_1 = ckx;

            sum += amp[k]*skx; 
        }

        buff[n] = sum;
        cx_1 = cx;
        sx_1 = sx;
        cx = cx_1*cwp - sx_1*swp;
        sx = sx_1*cwp + cx_1*swp;
    }
}

/*
    m5 is a combination of m2 and m4_rs

    sin((k+1)x) = 2sin(kx)cos(x) - sin((k-1)x
*/
static void
m5(float *amp, int K, float *buff, int N)
{
    float p, t, f, w, swp, cwp, st, ct, ckx, ckx_1, skx, skx_1, skx_2, cx, cx_1, sx, sx_1, sum, cx2;

    p = 1.0f/48000.0f;
    t = 0.0f;
    f = 1.0f/(2.0f*PI);
    w = 2 * PI * f;
    
    swp = sinf(w*p);
    cwp = cosf(w*p);
    st = sinf(t);
    ct = cosf(t);

    // TODO: change this to work when frequency changes
    // To change frequency, simply update cwp and swp in the loop
    // cwp = cosf(w2), swp = sinf(w2)
    // TODO: suppot changing amp[k] on a sample-by-sample basis
    // Common case:
    // Changing f often
    // Changing amp[k] not so much.

    /*
        sin(-x) = -sin(x)

        k
        -2  sin(-2x) = -sin(2x) = -(2sin(x)cos(x) - sin(0)) = -(2sin(x)cos(x))
        -1  sin(-x) = -sin(x)
         0  sin(0x) = 0.0f
         1  sin(x) = sx 

        computing 2*cos(x)    
    */
    cx = ct;
    sx = st;
    for (int n = 0; n < N; n++) {
        skx = 0.0f;
        skx_1 = -sx;
        skx_2 = -2.0f*cx*sx;
        sum = 0.0f;

        cx2 = 2.0f*cx;

        for (int k = 0; k < K; k++) {
            sum += amp[k]*skx; 

            skx_2 = skx_1;
            skx_1 = skx;
            skx = cx2*skx_1 - skx_2;
        }

        buff[n] = sum;

        cx_1 = cx;
        sx_1 = sx;
        cx = cx_1*cwp - sx_1*swp;
        sx = sx_1*cwp + cx_1*swp;
    }
}

static void
m4_r_sse_opt(float *amp, int K, float *buff, int N)
{
    float p = 1.0f/48000.0f;
    float t = 0.0f;
    float f = 1.0f/(2.0f*PI);
    float w = 2 * PI * f;
    
    float s4wp = sinf(4.0f*w*p);
    float c4wp = cosf(4.0f*w*p);
    float s4t = sinf(4.0f*t);
    float c4t = cosf(4.0f*t);
    float c4x, s4x, c4x_1, s4x_1;

    __m128 s_amp[K/4];
/*
    for (int k4 = 0; k4 < K/4; k4++) {
        s_amp[k4] = _mm_loadu_ps(amp+4*k4);
    }
*/

    __m128 s_skx_1;
    __m128 s_ckx_1;
    __m128 s_skx;
    __m128 s_ckx;
    __m128 s_sum;
    __m128 out;
    __m128 s_s4x;
    __m128 s_c4x;

    c4x_1 = c4t;
    s4x_1 = s4t;
    c4x = c4t;
    s4x = s4t;

    for (int n = 0; n < N; n++) {
        float x = w*(t+n*p);

        s_s4x = _mm_set1_ps(s4x);
        s_c4x = _mm_set1_ps(c4x);

        s_skx = _mm_setr_ps(0.0f, sinf(x), sinf(2.0f*x), sinf(3.0f*x));
        s_ckx = _mm_setr_ps(1.0f, cosf(x), cosf(2.0f*x), cosf(3.0f*x));

        s_sum = _mm_setzero_ps();
        for (int k4 = 0; k4 < K/4; k4++) {
            s_sum = _mm_add_ps(s_sum, _mm_mul_ps(s_amp[k4], s_skx));

            s_skx_1 = s_skx;
            s_ckx_1 = s_ckx; 
            s_skx = _mm_add_ps(_mm_mul_ps(s_skx_1, s_c4x), _mm_mul_ps(s_ckx_1, s_s4x));
            s_ckx = _mm_sub_ps(_mm_mul_ps(s_ckx_1, s_c4x), _mm_mul_ps(s_skx_1, s_s4x));
        }
        out = s_sum;
        out = _mm_hadd_ps(out, out);
        out = _mm_hadd_ps(out, out);
        _mm_store_ss(buff+n, out);

        c4x_1 = c4x;
        s4x_1 = s4x;
    
        c4x = c4x_1*c4wp - s4x_1*s4wp;
        s4x = s4x_1*c4wp + c4x_1*s4wp;
    }
}

static void
m4_r_sse(float *amp, int K, float *buff, int N)
{
    float p = 1.0f/48000.0f;
    float t = 0.0f;
    float f = 1.0f/(2.0f*PI);
    float w = 2 * PI * f;
    
    float s4wp = sinf(4.0f*w*p);
    float c4wp = cosf(4.0f*w*p);
    float s4t = sinf(4.0f*t);
    float c4t = cosf(4.0f*t);
    float ckx[4], ckx_1[4];
    float skx[4], skx_1[4];
    float c4x[N];
    float s4x[N];

    __m128 s_amp[K/4];

    for (int k4 = 0; k4 < K/4; k4++) {
        s_amp[k4] = _mm_loadu_ps(amp+4*k4);
    }


    // TODO: change this to work when frequency changes
    // To change frequency, simply update cwp and swp in the loop
    // cwp = cosf(w2), swp = sinf(w2)
    // TODO: suppot changing amp[k] on a sample-by-sample basis
    // Common case:
    // Changing f often
    // Changing amp[k] not so much.

    c4x[0] = c4t;
    s4x[0] = s4t;
    for (int n = 1; n < N; n++) {
        // float x = w*(t+n*p);
        // 4x = 4*w*(t+n*p)
        // 4*w*(t+(n+1)*p) = 4w(t+np+p) = 4w(t+np)+4wp
        c4x[n] = c4x[n-1]*c4wp - s4x[n-1]*s4wp;
        s4x[n] = s4x[n-1]*c4wp + c4x[n-1]*s4wp;
    }

    __m128 s_skx_1;
    __m128 s_ckx_1;
    __m128 s_skx;
    __m128 s_ckx;
    __m128 s_sum;

    for (int n = 0; n < N; n++) {
        skx_1[0] = 0.0f;
        ckx_1[0] = 1.0f;
        float x = w*(t+n*p);


        float sum[4] = {0};
        for (int i = 0; i < 4; i++) {
            skx_1[i] = sinf(i*x);
            ckx_1[i] = cosf(i*x);
            sum[i] += amp[i]*skx_1[i];
        }

        /* For some reason, these aren't the same? */
        s_skx_1 = _mm_setr_ps(0.0f, sinf(x), sinf(2.0f*x), sinf(3.0f*x));
        s_ckx_1 = _mm_setr_ps(1.0f, cosf(x), cosf(2.0f*x), cosf(3.0f*x));

        s_sum = _mm_mul_ps(s_amp[0], s_skx_1);
        __m128 s_s4x = _mm_set1_ps(sinf(4.0f*x));
        __m128 s_c4x = _mm_set1_ps(cosf(4.0f*x));

        float a0_skx[4];
        _mm_storeu_ps(a0_skx, s_skx);
        //printf("\n\nn: %d k4: %d\n", n, 0);
        for (int i = 0; i < 4; i++) {
            //printf("%f : s %f\n", skx[i], a0_skx[i]);
        }

        for (int k4 = 1; k4 < K/4; k4++) {
            for (int i = 0; i < 4; i++) {
                int k = k4*4 + i;
                // sin((k+4)x) = sin(kx+4x) instead of sin((k+1)x)
                skx[i] = skx_1[i]*c4x[n] + ckx_1[i]*s4x[n];
                ckx[i] = ckx_1[i]*c4x[n] - skx_1[i]*s4x[n];
                skx_1[i] = skx[i];
                ckx_1[i] = ckx[i];
                sum[i] += amp[k]*skx[i]; 
            }


            s_skx = _mm_add_ps(_mm_mul_ps(s_skx_1, s_c4x), _mm_mul_ps(s_ckx_1, s_s4x));
            s_ckx = _mm_sub_ps(_mm_mul_ps(s_ckx_1, s_c4x), _mm_mul_ps(s_skx_1, s_s4x));

            float a_skx[4];
            _mm_store_ps(a_skx, s_skx);
            // printf("\n\nn: %d k4: %d\n", n, k4);
            for (int i = 0; i < 4; i++) {
                // printf("%f : s %f\n", skx[i], a_skx[i]);
            } 

            // s_skx = _mm_setr_ps(sinf(4*k4*x), sinf((4*k4+1)*x), sinf((4*k4+2)*x), sinf((4*k4+3)*x));
            
            s_skx_1 = s_skx;
            s_ckx_1 = s_ckx; 
            s_sum = _mm_add_ps(s_sum, _mm_mul_ps(s_amp[k4], s_skx));
        }

        __m128 out = s_sum;
        out = _mm_hadd_ps(out, out);
        out = _mm_hadd_ps(out, out);
        _mm_store_ss(buff+n, out);
    }
}

/*
    Problem in parallelizng:
    Each n relies on a previous n
    Each k relies on previous k

    Sol 1:
    Compute all n values
    Parallel compute k for different n's

    Sol 2:
    Compute all k values
    Parallel compute n for different k's
*/
static void
m2_par(float *amp, int K, float *buff, int N)
{
    float p = 1.0f/48000.0f;
    float sin_p = sinf(p);
    float cos_p = cosf(p);

    float cos_t[N];
    float sin_t[N];

    __m128 sse_cos_t[N/4];
    __m128 sse_sin_t[N/4];

    float sin_k[K];

    /* Calculate Initial 4 values for vector. */
    float sse_sin_init[4];
    float sse_cos_init[4];

    sse_sin_init[0] = 0.0f;
    sse_sin_init[1] = sinf(p);

    sse_cos_init[0] = 1.0f;
    sse_cos_init[1] = cosf(p);

    float cos_p2 = 2.0f*sse_cos_init[1];

    sse_sin_init[2] = sse_sin_init[1]*cos_p2;
    sse_sin_init[3] = sse_sin_init[2]*cos_p2 - sse_sin_init[1];

    sse_cos_init[2] = cos_p2*sse_cos_init[1] - 1.0f;
    sse_cos_init[3] = cos_p2*sse_cos_init[2] - sse_cos_init[1];

    float s_cos_4p = cos_p2*sse_cos_init[3] - sse_cos_init[2];
    float s_sin_4p = cos_p2*sse_sin_init[3] - sse_sin_init[2];

    /* Load Initial Vectors */
    sse_sin_t[0] = _mm_loadu_ps(sse_sin_init);
    sse_cos_t[0] = _mm_loadu_ps(sse_cos_init);

    /* Load multiplication factors */
    __m128 cos_4p = _mm_load1_ps(&s_cos_4p);
    __m128 sin_4p = _mm_load1_ps(&s_sin_4p);

    for (int k = 0; k < K; k++) {
        /* TODO: Parallelize */
        /*
            Parallel version:

            Assuming 128bit vectors
            First vector is seeded with values of:
            v[0] = {sin(0p) sin(1p) sin(2p) sin(3p)}
            v[1] = {sin(4p) sin(5p) sin(6p) sin(7p)}
            Each vector we increase the angle by 4p
            We have a vector for sin and one for cos

            For each vector we must compute: 
            sin(np + 4p) = sin(np)cos(4p) + sin(4p)cos(np)
            cos(np + 4p) = cos(np)cos(4p) - sin(np)sin(4p)

            requires 2 mul, 1 add or sub
        */
        /* Init first vector */
        // sse_sin_t[];
        // sse_cos_t[];

        for (int n = 1; n < N/4; n++) {
            // mul_ps or mul_ss?
            /* Compute sin(np)*cos(4p) */
            __m128 snp_c4p = _mm_mul_ps(sse_sin_t[n-1], cos_4p); 

            /* Compute sin(4p)*cos(np) */
            __m128 s4p_cnp = _mm_mul_ps(sse_cos_t[n-1], sin_4p);

            /* Compute cos(4p)*cos(np) */
            __m128 c4p_cnp = _mm_mul_ps(sse_cos_t[n-1], cos_4p);

            /* Compute sin(4p)*sin(np) */
            __m128 s4p_snp = _mm_mul_ps(sse_sin_t[n-1], sin_4p);

            /* Compute sin(np)*cos(4p) + sin(4p)*cos(np) */
            sse_sin_t[n] = _mm_add_ps(snp_c4p, s4p_cnp);

            /* Compute cos(np)*cos(4p) - sin(4p)*sin(np) */
            sse_cos_t[n] = _mm_sub_ps(c4p_cnp, s4p_snp);
        }
        for (int n = 0; n < N/4; n++) {
            _mm_storeu_ps(buff+n*4, sse_sin_t[n]);
        }
    }

#if 0
    /* TODO: Parallelize */
    float c2;
    for (int n = 0; n < N; n++) {
        sin_k[0] = 0.0f;
        sin_k[1] = sin_t[n];
        buff[n] = amp[1]*sin_t[n]; // sin_k[0] + sin_k[1]
        c2 = 2.0f*cos_t[n];

        /* TODO: Parallelize */
        for (int k = 2; k < K; k++) { 
            sin_k[k] = c2*sin_k[k-1] - sin_k[k-2];
            buff[n] += amp[k]*sin_k[k];
        }
    }
#endif
}

/*
    Read this documentation, the method is not very understandable without it.

    Identities:
    sin(a+b) = sin(a)cos(b) + sin(b)cos(a)
    cos(a+b) = cos(a)cos(b) - sin(a)sin(b)
    s//in(nx) = 2sin((n-1)x)cos(x) - sin((n-2)x)
    cos(nx) = 2cos((n-1)x)cos(x) - cos((n-2)x)

    fs = sample rate
    n = sample number
    k = harmonic
    t = starting time of buffer

    Tables:
    S1[n] = sin(n/fs)
    C1[n] = cos(n/fs)
    S2[n] = sin(t + n/fs)
    C2[n] = cos(t + n/fs)
    S3[k,n] = sin(k(t + n/fs))
    C3[k,n] = cos(k(t + n/fs))

    Identity Recurrences:
    sin(n/fs) = 2sin((n-1)/fs)cos(1/fs) - sin((n-2)/fs)
    cos(n/fs) = 2cos((n-1)/fs)cos(1/fs) - cos((n-2)/fs)
    sin(t+n/fs) = sin(t)cos(n/fs) + sin(n/fs)cos(t)
    cos(t+n/fs) = cos(t)cos(n/fs) - sin(t)sin(n/fs)
    sin(k(t+n/fs)) = 2sin((k-1)(t+n/fs))cos(t+n/fs) - sin((k-2)(t+n/fs))
    cos(k(t+n/fs)) = 2cos((k-1)(t+n/fs))cos(t+n/fs) - cos((k-2)(t+n/fs))

    Table Recurrences:
    S1[n] = 2*S1[n-1]*C1[1] - S1[n-2]
    C1[n] = 2*C1[n-1]*C1[1] - C1[n-2]
    S1[0] = 0
    C1[0] = 1
    S1[1] = sin(1/fs)
    C1[1] = cos(1/fs)

    S2[n] = S2[0]*C1[n] + S1[n]*C2[0]
    C2[n] = C2[0]*C1[n] - S2[0]*S1[n]
    S2[0] = sin(t)
    C2[0] = cos(t)

    S3[k,n] = 2*S3[k-1,n]*C2[n] - S3[k-2][n]
    C3[k,n] = 2*C3[k-1,n]*C2[n] - S3[k-2][n]
    S3[0,n] = 0
    C3[0,n] = 1
    S3[1,n] = S2[n]
    C3[1,n] = C2[n]
    S3[0,0] = 0
    C3[0,0] = 1

    O(NK)
    Actually, this is the same method as m2?
*/

static void
m3(float *amp, int K, float *buff, int N)
{
    float fs = 48000.0f;
    float s1_n, c1_n, s1_0, c1_0, s1_n1, c1_n1, s1_n2, c1_n2, s1_1, c1_1;
    float s2_n, s2_0, c2_n, c2_0;
    float s3_k, c3_k, s3_k1, s3_k2, s3_0, s3_1, c3_0, c3_1;

    float n_s3_k;

    s1_0 = 0.0f;
    c1_0 = 1.0f;
    s1_1 = sinf(1.0/fs);
    c1_1 = sinf(1.0/fs);

    s2_0 = 0.0f; //sin(t), t=0 in this test
    c2_0 = 1.0f; //cos(t), t=0 in this test 

    s3_0 = 0.0f;
    c3_0 = 1.0f;

    c1_n1 = c1_1;
    c1_n2 = c1_0;
    s1_n1 = s1_1;
    s1_n2 = s1_0;

    for (int n = 0; n < N; n++) {
        c1_n = 2.0f*c1_n1*c1_1 - c1_n2;
        s1_n = 2.0f*s1_n1*c1_1 - s1_n2;

        c2_n = c2_0*c1_n - s2_0*s1_n;
        s2_n = s2_0*c1_n + s1_n*c2_0;

        buff[n] = 0.0f;
       
        /* Correct? */ 
        s3_k1 = s2_n;
        s3_k2 = 0.0f;
        for (int k = 0; k < K; k++) {
            s3_k = 2.0f*s3_k1*c2_n - s3_k2;

            buff[n] += amp[k]*s3_k;
            s3_k2 = s3_k1;
            s3_k1 = s3_k;
        }

        c1_n1 = c1_n;
        c1_n2 = c1_n1;
        s1_n1 = s1_n;
        s1_n2 = s1_n1;
    }
}

int
main(int argc, char *argv[])
{

    /* Flush to zero, as denormalized floats can cause big slowdowns. */
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

    if (argc != 2) {
        printf("Usage: %s iterations\n", argv[0]);
        return 1;
    }

    /* Wavetable setup not included in benchmark, as it would be run only on sample rate change. */
    G_SINETABLE_LEN = 1024;
    G_SINETABLE = malloc(sizeof(float)*G_SINETABLE_LEN);
    for (int i = 0; i < G_SINETABLE_LEN; i++) {
        G_SINETABLE[i] = sinf(2*PI*i/(double)G_SINETABLE_LEN);
    }

    int iters = atoi(argv[1]);
    bmark(direct, "direct", iters);
    bmark(m2, "m2", iters); 
    bmark(m2_par, "m2_par", iters); 
    bmark(m2_par2, "m2_par2", iters); 
    bmark(m3, "m3", iters); 
    bmark(m4, "m4", iters); 
    bmark(m4_arr, "m4_arr", iters); 
    bmark(m4_sse, "m4_sse", iters); 
    bmark(m4_r, "m4_r", iters); 
    bmark(m4_rs, "m4_rs", iters); 
    bmark(m4_r_sse, "m4_r_sse", iters); 
    bmark(m4_r_sse_opt, "m4_r_sse_opt", iters); 
    bmark(m4_sse_opt, "m4_sse_opt", iters); 
    bmark(m4_sse_opt2, "m4_sse_opt2", iters); 
    bmark(wave_table, "wave_table", iters); 
    bmark(m5, "m5", iters); 

    free(G_SINETABLE);
    return 0;
}
