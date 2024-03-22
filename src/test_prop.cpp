#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
#include <ctime>
//#include <emmintrin.h>

template <class T>
T findMax(const std::vector<T>& data){
    T __curMax = 0;
    for(unsigned int i = 0; i < data.size(); ++i){
        if(data[i] > __curMax){
            __curMax = data[i];
        }
    }
    return __curMax;
}

template <class T>
std::vector<T> buildReference(const IntegerMatrix& reference, int shift = 2){
    int ncolR = reference.ncol();
    int nrowR = reference.nrow();
    
    // Build integer representations
    std::vector<T> __refValues(nrowR, 0);
    for(int i = 0; i < nrowR; ++i){
        // Binary magic
        for(int j = 0; j < ncolR; ++j){
            __refValues[i] <<= shift;        
            __refValues[i] ^= reference(i,j);
        }
    }
    
     // Find max
    int __curMax = findMax<int>(__refValues);
    
    // Build lookup map
    std::vector<int> __ref(__curMax+1);
     for(int i = 0; i < nrowR; ++i){
         __ref[__refValues[i]] = i;
     }
     
    return __ref;
}

// #define ROUND_DOWN(x, s) ((x) & ~((s)-1))
//  // commenting out threads = omp_get_num_threads();
// // [[Rcpp::export]]
// NumericMatrix combinationRefs(const IntegerMatrix& x, const IntegerMatrix& cols, const IntegerMatrix& reference, const NumericVector& prop, int shiftSize = 2, int threads = -1, bool verbose = false){
//     const int nrow = x.nrow();
//     // const int ncol = x.ncol();
//     const int nrow2 = cols.nrow();
//     const int ncol2 = cols.ncol();
//     const int ncolR = reference.ncol();
//     const int nrowR = reference.nrow();
    
//     if(threads <= 0)
//         //threads = omp_get_num_threads();
// 	threads=1;
    
//     // Prepare reference
//     if(verbose)
//         Rcpp::Rcout << "Preparing reference..." << std::endl;
    
//     // Build ref
//     std::vector<int> __ref = buildReference<int>(reference, shiftSize);
    
//     // Preallocate memory for matrix
//     if(verbose)
//         Rcpp::Rcout << "Allocate matrix (" << ncol2 << ", " << nrowR << ": " << (unsigned long long)ncol2*nrowR << " entries)" << std::endl;
    
//     NumericMatrix ret(ncol2, nrowR);
   
//    if(verbose)
//         Rcpp::Rcout << "Begin (" << threads << " threads)..." << std::endl;

//     std::clock_t start;
//     start = std::clock();
    
//     // Can't increase cache-line efficiency because of how R structures are designed
//     // Runtime O(n*(mr+s)) ~~~ O(n(r+s))
//     //pragma omp parallel for num_threads(threads)
//     for(int k = 0; k < ncol2; ++k){
//        int* vals = 0;
//        if(posix_memalign((void**)&vals, 16, nrow * sizeof(int)))
// 	 return 0;
//        memset(&vals[0], 0, nrow * sizeof vals[0]); // Stupid OpenMP limitation
        
//         // int positions[nrow2+1];
// 	int *positions = new int[nrow2+1];
//   // char filename1char[filenameLength];
//   // char *filename1 = new char[filenameLength];
//  for(int i = 0; i < nrow2; ++i)
//             positions[i] = cols(i,k)-1;

//          // Frequency
//          for(int j = 0; j < ncolR; ++j){
//             const int* p = &x[positions[j]*nrow]; // Stupid R limitation
//             const int stepsize = 4;
//             int i = 0;
//             // Binary magic using SSE vectorization and manual loop unrolling
//             // Access data linearly
//             for(; i < ROUND_DOWN(nrow, stepsize); i+=stepsize){
//                 __m128i v0 = _mm_load_si128((__m128i*)&vals[i]);
//                 __m128i P = _mm_loadu_si128((__m128i*)&p[i]); // Stupid R limitation
//                 v0 = _mm_slli_epi32(v0, shiftSize);
//                 v0 = _mm_xor_si128 (v0, P);
//                 _mm_store_si128((__m128i*)&vals[i], v0);
//             }
//             for(; i < nrow; i++){
//                 vals[i] <<= shiftSize;
//                 vals[i] ^= p[i];
//             }
//         }
            
//             std::vector<double> sums(nrowR, 0); // Stupid OpenMP limitation
//          for(int i = 0; i < nrow; ++i)
//             sums[__ref[vals[i]]] += prop[i]; // Chain dependency, cannot unroll
        
//         // Fill probabilities into output matrix
//         for(int i = 0; i < nrowR; ++i)
//             ret(k, i) = sums[i];
        
//         // Cleanup
//         free(vals);
//     }
    
//     if(verbose){
//         Rcpp::Rcout << "Finished " << (unsigned long long)ncol2*(nrow*ncolR+nrow+nrowR) << " operations in body" << std::endl;
//         Rcpp::Rcout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
//     }
    
//     return(ret);
// }

// // [[Rcpp::export]]
// IntegerMatrix combination2(IntegerMatrix x, IntegerMatrix cols) {
//     int nrow = x.nrow();
//     int ncol = x.ncol();
//     int nrow2 = cols.nrow();
//     int ncol2 = cols.ncol();
//     int val = 0;
    
//     IntegerMatrix ret(nrow,ncol2);
//     std::vector<int> positions(3,0);
     
//     for(int k = 0; k < ncol2; ++k){
// 		for(int i = 0; i < nrow2; ++i)
//             positions[i] = cols(i,k)-1;
//         for(int i = 0; i < nrow; ++i){ // Rows
//             val = 0;
//             for(int j = 0; j < 3; ++j){ // Selected columns
//                 val <<= 2;
//                 val ^= x(i,positions[j]) << 0;
//             }
//             ret(i,k) = val;
//         }
//     }
//     return(ret);
// }
