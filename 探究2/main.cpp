//g++ avx_pi.cpp -mavx -O2
#include <iostream>
#include <ctime>
#include <cmath>
#include <x86intrin.h>

double piLeibnizAVX(size_t dt);

double piLeibnizNaive(size_t dt);

double piFourierAVX(size_t dt);

double piFourierNaive(size_t dt);

double Ramanujan(int iter);

using namespace std;

//正常的逐个累加运算
double piIntegrationNaive(size_t dt) {
    double pi = 0.0;
    double delta = 1.0 / dt;
    for (size_t i = 0; i < dt; i++) {
        double x = (double) i / dt;
        pi += delta / (1 + x * x);
    }
    return pi * 4.0;
}

//利用avx256指令集
double piIntegrationAVX(size_t dt) {
    double pi = 0.0;
    double delta = 1.0 / dt;
    __m256d sgn_p_vec, delta_vec, base_vec, temp_vec, res_vec;
    sgn_p_vec = _mm256_set1_pd(1.0);
    delta_vec = _mm256_set1_pd(delta);
    base_vec = _mm256_set_pd(delta * 3, delta * 2, delta, 0.0);
    res_vec = _mm256_setzero_pd();
    for (size_t i = 0; i < dt - 4; i += 4) {
        temp_vec = _mm256_set1_pd(i * delta);
        temp_vec = _mm256_add_pd(temp_vec, base_vec);
        temp_vec = _mm256_mul_pd(temp_vec, temp_vec);
        temp_vec = _mm256_add_pd(sgn_p_vec, temp_vec);
        temp_vec = _mm256_div_pd(delta_vec, temp_vec);
        res_vec = _mm256_add_pd(res_vec, temp_vec);
    }
    double res[4] __attribute__((aligned(32)));
    _mm256_store_pd(res, res_vec);
    pi += res[0] + res[1] + res[2] + res[3];
    return pi * 4.0;
}

int test(long long ratio);

int main() {
    long long ratio = 1;
    for (int i = 0; i < 1; i++) {
        test(ratio);
        ratio *= 10;
    }
}


int test(long long ratio) {
    clock_t start, end;
    size_t dt = 1e10 / ratio;
    double result1, result2;

    long long time_sum = 0;
    for (int i = 0; i < 100000000; i++) {
        start = clock();
        result1 = Ramanujan(1);
        end = clock();
        time_sum += end - start;
    }
    cout << "Ramanujan:\n" << result1 << endl << time_sum << endl << endl;

   //普通函数计时
   size_t dt1 = dt/2;
   for (int i = 0; i < 1 * ratio; i++) {
       start = clock();
       result1 = piFourierNaive(dt);
       end = clock();
       time_sum += end - start;
   }
   cout << "naive:\n" << result1 << endl << time_sum << endl << endl;
//
//    //avx256计时
//    time_sum = 0;
//    for (int i = 0; i < 1 * ratio; i++) {
//        start = clock();
//        result2 = piFourierAVX(dt);
//        end = clock();
//        time_sum += end - start;
//    }
//    cout << "avx256:\n" << result2 << endl << time_sum << endl << endl;


    return 0;
}

double Ramanujan(int iter) {
    long long factorial[16];
    factorial[0] = 1;
    factorial[1] = 1;
    for (int i = 2; i <= iter * 4; i++) {
        factorial[i] = factorial[i - 1] * i;
    }

    double res = 0, temp = 0;
    for (int i = 0; i < iter; i++) {
        temp = factorial[4 * i];
        temp /= std::pow(factorial[i], 4);
        temp /= std::pow(396, 4 * i);
        temp *= 26390 * i + 1103;
        res += temp;
    }
    res *= 2 * std::sqrt(2);
    res /= 99 * 99;
    return 1.0 / res;
}

double piFourierNaive(size_t dt) {
    double res = 0;
    for (size_t i = 1; i < dt; i += 2) {
        res += 1.0 / (i * i);
    }

    return std::sqrt(res * 8);
}

double piFourierAVX(size_t dt) {
    __m256d base, temp, res_vec, sgn_p;
    sgn_p = _mm256_set1_pd(1);
    base = _mm256_set_pd(0, 2, 4, 6);
    res_vec = _mm256_set1_pd(0);
    double res = 0;
    for (size_t i = 1; i < dt; i += 8) {
        temp = _mm256_set1_pd(i);
        temp = _mm256_add_pd(base, temp);
        temp = _mm256_mul_pd(temp, temp);
        temp = _mm256_div_pd(sgn_p, temp);
        res_vec = _mm256_add_pd(res_vec, temp);
    }
    double tmp[4] __attribute__((aligned(32)));
    _mm256_store_pd(tmp, res_vec);
    res += tmp[0] + tmp[1] + tmp[2] + tmp[3];
    return std::sqrt(res * 8);
}

double piLeibnizNaive(size_t dt) {
    double res = 0;
    for (size_t i = 1; i < dt; i += 4) {
        res += 1.0 / i;
    }
    for (size_t i = 3; i < dt; i += 4) {
        res -= 1.0 / i;
    }

    return res * 4;
}

double piLeibnizAVX(size_t dt) {
    __m256d base, sgn_p, sgn_n, temp, res_vec;
    base = _mm256_set_pd(0, 4, 8, 12);
    sgn_p = _mm256_set1_pd(1.0);
    sgn_n = _mm256_set1_pd(-1.0);
    res_vec = _mm256_set1_pd(0);
    double res = 0;
    for (size_t i = 1; i < dt; i += 16) {
        temp = _mm256_set1_pd(i);
        temp = _mm256_add_pd(base, temp);
        temp = _mm256_div_pd(sgn_p, temp);
        res_vec = _mm256_add_pd(res_vec, temp);
    }
    double tmp[4] __attribute__((aligned(32)));
    _mm256_store_pd(tmp, res_vec);
    res += tmp[0] + tmp[1] + tmp[2] + tmp[3];
    res_vec = _mm256_set1_pd(0);

    for (size_t i = 3; i < dt; i += 16) {
        temp = _mm256_set1_pd(i);
        temp = _mm256_add_pd(base, temp);
        temp = _mm256_div_pd(sgn_n, temp);
        res_vec = _mm256_add_pd(res_vec, temp);
    }
    _mm256_store_pd(tmp, res_vec);
    res += tmp[0] + tmp[1] + tmp[2] + tmp[3];
    return res * 4;
}
