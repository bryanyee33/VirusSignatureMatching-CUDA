#include "kseq/kseq.h"
#include "common.h"

#define SAMPLE_MAX_LEN 200000
#define SIGNATURE_MAX_LEN 10000

// SAMPLE_SMEM_SIZE = SAMPLE_MAX_LEN / 2 + SIGNATURE_MAX_LEN
#define SAMPLE_SMEM_SIZE 110000

void check_cuda_errors()
{
    cudaError_t rc;
    rc = cudaGetLastError();
    if (rc != cudaSuccess)
    {
        printf("Last CUDA error %s\n", cudaGetErrorString(rc));
    }

}

__global__ void getMatches(const char* __restrict d_sample_seq, const char* __restrict d_sample_qual, const char* __restrict d_signature_seq,
        double* __restrict match_scores, unsigned short* __restrict match_idx, int* __restrict match_count) {

    // Half of sample so it will always fit in shared memory; A size of SIGNATURE_MAX_LEN needs to be copied twice
    extern __shared__ char s_sample[];
    const char *sample = &d_sample_seq[SAMPLE_MAX_LEN * blockIdx.x];
    const char *signature = &d_signature_seq[SIGNATURE_MAX_LEN * threadIdx.x];

    // Copy 1st half + SIGNATURE_MAX_LEN of sample into shared memory
    int num_per_thread = (SAMPLE_SMEM_SIZE + blockDim.x - 1) / blockDim.x; // Ceiling
    if (threadIdx.x == blockDim.x - 1) { // Last thread
        memcpy(&s_sample[threadIdx.x * num_per_thread], &sample[threadIdx.x * num_per_thread],
               SAMPLE_SMEM_SIZE - num_per_thread * threadIdx.x);
    } else {
        memcpy(&s_sample[threadIdx.x * num_per_thread], &sample[threadIdx.x * num_per_thread],
               num_per_thread);
    }

    // if (threadIdx.x == 0) {
    //     memcpy(s_sample, sample, SAMPLE_SMEM_SIZE);
    // }
     __syncthreads();

    for (int i = 0; i < (SAMPLE_MAX_LEN >> 1); ++i) {
        if (!s_sample[i]) { // Sample ended; no match
            return;
        }

        for (int j = 0; j < SIGNATURE_MAX_LEN; ++j) {
            // TODO: FIX WHEN SIGNATURE == MAX LEN
            if (!signature[j]) { // Signature ended; found match
                const char *qual = &d_sample_qual[SAMPLE_MAX_LEN * blockIdx.x + i];
                int tot = 0;
                for (int s = 0; s < j; ++s) {
                    tot += qual[s] - 33;
                }

                int idx = atomicAdd(match_count, 1); // Each idx will be unique (no need synchronisation)
                match_scores[idx] = (double)tot / j;
                match_idx[idx << 1] = blockIdx.x; // Sample idx
                match_idx[(idx << 1) + 1] = threadIdx.x; // Signature idx
                return;

            } else if (!s_sample[i + j] || s_sample[i + j] != 'N' && signature[j] != 'N' &&
                    s_sample[i + j] != signature[j]) { // No match
                break;
            }
        }
    }

    // Copy 2nd half of sample into shared memory
    num_per_thread = ((SAMPLE_MAX_LEN >> 1) + blockDim.x - 1) / blockDim.x; // Ceiling
    if (threadIdx.x == blockDim.x - 1) { // Last thread
        memcpy(&s_sample[threadIdx.x * num_per_thread], &sample[(SAMPLE_MAX_LEN >> 1) + threadIdx.x * num_per_thread],
               (SAMPLE_MAX_LEN >> 1) - num_per_thread * threadIdx.x);
        s_sample[(SAMPLE_MAX_LEN >> 1)] = 0; // Set next char to '\0' to indicate end
    } else {
        memcpy(&s_sample[threadIdx.x * num_per_thread], &sample[(SAMPLE_MAX_LEN >> 1) + threadIdx.x * num_per_thread],
               num_per_thread);
    }

    __syncthreads();

    for (int i = 0; i < (SAMPLE_MAX_LEN >> 1); ++i) {
        if (!s_sample[i]) { // Sample ended; no match
            return;
        }

        for (int j = 0; j < SIGNATURE_MAX_LEN; ++j) {
            // TODO: FIX WHEN SIGNATURE == MAX LEN
            if (!signature[j]) { // Signature ended; found match
                const char *qual = &d_sample_qual[SAMPLE_MAX_LEN * blockIdx.x + i];
                int tot = 0;
                for (int s = 0; s < j; ++s) {
                    tot += qual[s] - 33;
                }

                int idx = atomicAdd(match_count, 1); // Each idx will be unique (no need synchronisation)
                match_scores[idx] = (double)tot / j;
                match_idx[idx << 1] = blockIdx.x; // Sample idx
                match_idx[(idx << 1) + 1] = threadIdx.x; // Signature idx
                return;

            } else if (!s_sample[i + j] || s_sample[i + j] != 'N' && signature[j] != 'N' &&
                    s_sample[i + j] != signature[j]) { // No match
                break;
            }
        }
    }
}

void runMatcher(const std::vector<klibpp::KSeq>& samples, const std::vector<klibpp::KSeq>& signatures, std::vector<MatchResult>& matches) {
    int sample_arr_len = samples.size() * SAMPLE_MAX_LEN;
    int signature_arr_len = signatures.size() * SIGNATURE_MAX_LEN;

    // calloc() to initialise all to '\0'; set to max possible length for each type.
    char *h_sample_seq = (char*)calloc(sample_arr_len, sizeof(char));
    char *h_sample_qual = (char*)calloc(sample_arr_len, sizeof(char));
    char *h_signature_seq = (char*)calloc(signature_arr_len, sizeof(char));

    // Device variables
    char *d_sample_seq, *d_sample_qual, *d_signature_seq;
    int *match_count; // 1 int storing the number of matches

    // Pinned memory in host
    double *match_scores;      // [Score_1, Score2, ...]
    unsigned short *match_idx; // [Samp_idx_1, Sig_idx_1, Samp_idx_2, Sig_idx_2, ...]

    // Combine each string array into 1 for transfer
    for (int i = 0; i < samples.size(); ++i) {
        memcpy(&h_sample_seq[i * SAMPLE_MAX_LEN], &samples[i].seq[0], samples[i].seq.size());
        memcpy(&h_sample_qual[i * SAMPLE_MAX_LEN], &samples[i].qual[0], samples[i].qual.size());
    }
    for (int i = 0; i < signatures.size(); ++i) {
        memcpy(&h_signature_seq[i * SIGNATURE_MAX_LEN], &signatures[i].seq[0], signatures[i].seq.size());
    }

    // Use cudaMallocHost() to write straight to host, since only a small number of matches should be found
    cudaMallocHost(&match_scores, sizeof(double) * samples.size() * signatures.size()); // Max possible number of matches
    cudaMallocHost(&match_idx, sizeof(unsigned short) * 2 * samples.size() * signatures.size());

    // Test out zero-copy memory
    // cudaHostAlloc(&match_scores, sizeof(double) * samples.size() * signatures.size(), cudaHostAllocMapped); // Max possible number of matches
    // cudaHostAlloc(&match_idx, sizeof(unsigned short) * 2 * samples.size() * signatures.size(), cudaHostAllocMapped);
    cudaMalloc(&match_count, sizeof(int));
    cudaMemset(match_count, 0, 1); // Initialise match_count to 0

    cudaMalloc(&d_sample_seq, sample_arr_len);
    cudaMalloc(&d_sample_qual, sample_arr_len);
    cudaMalloc(&d_signature_seq, signature_arr_len);
    cudaMemcpy(d_sample_seq, h_sample_seq, sample_arr_len, cudaMemcpyHostToDevice);
    cudaMemcpy(d_sample_qual, h_sample_qual, sample_arr_len, cudaMemcpyHostToDevice);
    cudaMemcpy(d_signature_seq, h_signature_seq, signature_arr_len, cudaMemcpyHostToDevice);

    cudaFuncSetAttribute(getMatches, cudaFuncAttributeMaxDynamicSharedMemorySize, SAMPLE_SMEM_SIZE);
    // 1 Sample per block; Each thread in block corresponds to 1 signature
    getMatches<<<samples.size(), signatures.size(), SAMPLE_SMEM_SIZE>>>(d_sample_seq, d_sample_qual, d_signature_seq, match_scores, match_idx, match_count);
    check_cuda_errors();
    cudaDeviceSynchronize();

    int h_match_count;
    cudaMemcpy(&h_match_count, match_count, sizeof(int), cudaMemcpyDeviceToHost);
    for (int i = 0; i < h_match_count; ++i) {
        matches.emplace_back(MatchResult(samples[match_idx[i << 1]].name,
                             signatures[match_idx[(i << 1) + 1]].name,
                             match_scores[i]));
    }

    cudaFree(d_sample_seq);
    cudaFree(d_sample_qual);
    cudaFree(d_signature_seq);
    free(h_sample_seq);
    free(h_sample_qual);
    free(h_signature_seq);
}
