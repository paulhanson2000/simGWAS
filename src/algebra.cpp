#include <Rcpp.h>
using namespace Rcpp;
#include <ctime>

#define ROUND_DOWN(x, s) ((x) & ~((s)-1))

#ifdef __ARM_NEON
#include <arm_neon.h>
// [[Rcpp::export]]
NumericVector MatrixVector(const NumericMatrix& matrix, const NumericVector& vector, const bool verbose = false){
    const int nrow = matrix.nrow();
    const int ncol = matrix.ncol();
    NumericVector ret(nrow);
    const double* ref;

    std::clock_t start;
    start = std::clock();

    for(int i = 0; i < ncol; ++i){
        ref = &matrix[nrow*i];
        float64x2_t P = vdupq_n_f64(vector[i]);

        int j = 0;
        for(; j < ROUND_DOWN(nrow, 2); j+=2){ // Process two elements at a time
            float64x2_t v0 = vld1q_f64(&ref[j]); // Load 2 elements from ref into a NEON register
            float64x2_t r = vld1q_f64(&ret[j]); // Load 2 elements from ret into a NEON register
            v0 = vmulq_f64(v0, P); // Multiply elements by vector[i]
            v0 = vaddq_f64(v0, r); // Add elements to ret
            vst1q_f64(&ret[j], v0); // Store result back to ret
        }
        for(; j < nrow; ++j){ // Process remaining elements
            ret[j] += ref[j] * vector[i]; // Perform scalar multiplication
        }
    }

    if(verbose){
        Rcout << "Clocks/cycle: " << CLOCKS_PER_SEC << std::endl;
        const double timediff = (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
        Rcout << "Operations: " << ncol*nrow*2 << " --> " << (ncol*nrow*2)/timediff/1000 << " MFLOPS" << std::endl;
        Rcpp::Rcout << "Time: " << timediff << " ms" << std::endl;
    }
    return(ret);
}
#else // Intel SSE
#include <emmintrin.h>
// [[Rcpp::export]]
NumericVector MatrixVector(const NumericMatrix& matrix, const NumericVector& vector, const bool verbose = false){
    const int nrow = matrix.nrow();
    const int ncol = matrix.ncol();
    NumericVector ret(nrow);
    const double* ref;
    
    std::clock_t start;
    start = std::clock();
    
    for(int i = 0; i < ncol; ++i){
        ref = &matrix[nrow*i];
        // Load vector[i]
        __m128d P = _mm_set1_pd(vector[i]);
        
        int j = 0;
         for(; j < ROUND_DOWN(nrow, 4); j+=4){
            __m128d v0 = _mm_loadu_pd(&ref[j]); // Stupid R limitation
            __m128d v1 = _mm_loadu_pd(&ref[j+2]); // Stupid R limitation
            __m128d r = _mm_loadu_pd(&ret[j]); // Stupid R limitation
             __m128d r1 = _mm_loadu_pd(&ret[j+2]); // Stupid R limitation
            v0 = _mm_mul_pd(v0, P);
            v1 = _mm_mul_pd(v1, P);
            v0 = _mm_add_pd(v0, r);
            v1 = _mm_add_pd(v1, r1);
            _mm_storeu_pd(&ret[j], v0); // Stupid R limitation
             _mm_storeu_pd(&ret[j+2], v1); // Stupid R limitation
        }
        for(; j < nrow; ++j){
            ret[j] += ref[j]*vector[i];
        }
    }
    
    if(verbose){
        Rcout << "Clocks/cycle: " << CLOCKS_PER_SEC << std::endl;
        const double timediff = (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
        Rcout << "Operations: " << ncol*nrow*2 << " --> " << (ncol*nrow*2)/timediff/1000 << " MFLOPS" << std::endl;
         Rcpp::Rcout << "Time: " << timediff << " ms" << std::endl;
    }
    return(ret);
}
#endif
