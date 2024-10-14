#include "kseq/kseq.h"
#include "common.h"

#define SAMPLE_MAX_LEN 200000
#define SIGNATURE_MAX_LEN 10000

void check_cuda_errors() {
    cudaError_t rc = cudaGetLastError();
    if (rc != cudaSuccess) {
        printf("Last CUDA error %s\n", cudaGetErrorString(rc));
    }
}

__global__ void getMatches(const char* __restrict d_sample_seq, const char* __restrict d_sample_qual,
        const char* __restrict d_signature_seq, const char* __restrict d_signature_qual) {

    const char *sample = &d_sample_seq[SAMPLE_MAX_LEN * blockIdx.x];
    const char *signature = &d_signature_seq[SIGNATURE_MAX_LEN * threadIdx.x];

    for (int i = 0; i < SAMPLE_MAX_LEN; ++i) {
        for (int j = 0; j < SIGNATURE_MAX_LEN; ++j) {
            if (!signature[j]) { // Signature ended; found match
                //i to j-1
                
                printf("Sample: %d, Sig: %d, %c\n", blockIdx.x, threadIdx.x, signature[j - 1]);
                break;
            } else if (!sample[i + j] || sample[i + j] != 'N' && signature[j] != 'N' &&
                    sample[i + j] != signature[j]) { // No match
                break;
            }
        }
    }
    // const int id = threadIdx.x + blockIdx.x * blockDim.x;

    // __shared__ unsigned long long data[BLOCK_SIZE];

    // data[threadIdx.x] = input[id];
    // __syncthreads();

    // for (int i = 1; i < 256; i <<= 1) {
    //     if (threadIdx.x % (i << 1) == 0) {
    //         data[threadIdx.x] += data[threadIdx.x + i];
    //     }
    //     __syncthreads();
    // }

    // if (threadIdx.x == 0) {
    //     atomicAdd(&dResult, data[0]);
    // }
}

void runMatcher(const std::vector<klibpp::KSeq>& samples, const std::vector<klibpp::KSeq>& signatures, std::vector<MatchResult>& matches) {
    int sample_arr_len = samples.size() * SAMPLE_MAX_LEN;
    int signature_arr_len = signatures.size() * SIGNATURE_MAX_LEN;

    // calloc() to initialise all to '\0'; set to max possible length for each type.
    char *h_sample_seq = (char*)calloc(sample_arr_len, sizeof(char));
    char *h_sample_qual = (char*)calloc(sample_arr_len, sizeof(char));
    char *h_signature_seq = (char*)calloc(signature_arr_len, sizeof(char));
    char *h_signature_qual = (char*)calloc(signature_arr_len, sizeof(char));

    for (int i = 0; i < samples.size(); ++i) {
        memcpy(&h_sample_seq[i * SAMPLE_MAX_LEN], &samples[i].seq[0], samples[i].seq.size());
        memcpy(&h_sample_qual[i * SAMPLE_MAX_LEN], &samples[i].seq[0], samples[i].qual.size());
    }
    for (int i = 0; i < signatures.size(); ++i) {
        memcpy(&h_signature_seq[i * SIGNATURE_MAX_LEN], &signatures[i].seq[0], signatures[i].seq.size());
        memcpy(&h_signature_qual[i * SIGNATURE_MAX_LEN], &signatures[i].seq[0], signatures[i].qual.size());
    }
    
    char *d_sample_seq, *d_sample_qual, *d_signature_seq, *d_signature_qual;

    cudaMallocManaged(&d_sample_seq, sample_arr_len);
    cudaMallocManaged(&d_sample_qual, sample_arr_len);
    cudaMemcpy(d_sample_seq, h_sample_seq, sample_arr_len, cudaMemcpyHostToDevice);
    cudaMemcpy(d_sample_qual, h_sample_qual, sample_arr_len, cudaMemcpyHostToDevice);

    cudaMallocManaged(&d_signature_seq, signature_arr_len);
    cudaMallocManaged(&d_signature_qual, signature_arr_len);
    cudaMemcpy(d_signature_seq, h_signature_seq, signature_arr_len, cudaMemcpyHostToDevice);
    cudaMemcpy(d_signature_qual, h_signature_qual, signature_arr_len, cudaMemcpyHostToDevice);

    getMatches<<<samples.size(), signatures.size()>>>(d_sample_seq, d_sample_qual, d_signature_seq, d_signature_qual);
    check_cuda_errors();

    cudaFree(d_sample_seq);
    cudaFree(d_sample_qual);
    cudaFree(d_signature_seq);
    cudaFree(d_signature_qual);

    // unsigned long long result = 0.0f;
    // // Run several iterations to get an average measurement
    // for (unsigned long long i = 0; i < TIMING_ITERATIONS; i++)
    // {
    //     // Reset acummulated result to 0 in each run
    //     cudaMemcpyToSymbol(dResult, &result, sizeof(unsigned long long));
    //     func<<<gridDim, blockDim>>>(dValsPtr, N);
    // }

    // // cudaMemcpyFromSymbol will implicitly synchronize CPU and GPU
    // cudaMemcpyFromSymbol(&result, dResult, sizeof(unsigned long long));

    free(h_sample_seq);
    free(h_sample_qual);
    free(h_signature_seq);
    free(h_signature_qual);
}
