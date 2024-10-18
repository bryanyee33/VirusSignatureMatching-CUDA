#include "kseq/kseq.h"
#include "common.h"

#define SAMPLE_MAX_LEN 200000
#define SIGNATURE_MAX_LEN 10000

__global__ void getMatches(const char* __restrict d_sample_seq, const char* __restrict d_sample_qual, const char* __restrict d_signature_seq,
        double* __restrict match_scores, unsigned short* __restrict match_idx, int* __restrict match_count) {

    const char *sample = &d_sample_seq[SAMPLE_MAX_LEN * blockIdx.x];
    const char *signature = &d_signature_seq[SIGNATURE_MAX_LEN * threadIdx.x];

    for (int i = 0; i < SAMPLE_MAX_LEN; ++i) {
        if (!sample[i]) { // Sample ended; no match
            return;
        }

        for (int j = 0; j < SIGNATURE_MAX_LEN + 1; ++j) { // +1 iteration for when signature len == SIGNATURE_MAX_LEN
            if (j == SIGNATURE_MAX_LEN || !signature[j]) { // Signature ended; found match
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

            } else if (!sample[i + j] || sample[i + j] != 'N' && signature[j] != 'N' &&
                    sample[i + j] != signature[j]) { // No match
                break;
            }
        }
    }
}

void runMatcher(const std::vector<klibpp::KSeq>& samples, const std::vector<klibpp::KSeq>& signatures, std::vector<MatchResult>& matches) {
    int sample_arr_len = samples.size() * SAMPLE_MAX_LEN;
    int signature_arr_len = signatures.size() * SIGNATURE_MAX_LEN;

    // Device variables
    char *d_sample_seq, *d_sample_qual, *d_signature_seq;
    int *match_count; // 1 int storing the number of matches
    cudaMalloc(&d_sample_seq, sample_arr_len);
    cudaMemset(d_sample_seq, 0, sample_arr_len); // Initialise arrays to '\0'

    cudaMalloc(&d_sample_qual, sample_arr_len);
    cudaMemset(d_sample_qual, 0, sample_arr_len);

    cudaMalloc(&d_signature_seq, signature_arr_len);
    cudaMemset(d_signature_seq, 0, signature_arr_len);

    cudaMalloc(&match_count, sizeof(int));
    cudaMemset(match_count, 0, 1); // Initialise match_count to 0

    // Pinned memory
    double *match_scores;      // [Score_1, Score2, ...]
    unsigned short *match_idx; // [Samp_idx_1, Sig_idx_1, Samp_idx_2, Sig_idx_2, ...]

    for (int i = 0; i < samples.size(); ++i) {
        cudaMemcpy(&d_sample_seq[i * SAMPLE_MAX_LEN], samples[i].seq.c_str(), samples[i].seq.size(), cudaMemcpyHostToDevice);
        cudaMemcpy(&d_sample_qual[i * SAMPLE_MAX_LEN], samples[i].qual.c_str(), samples[i].qual.size(), cudaMemcpyHostToDevice);
    }
    for (int i = 0; i < signatures.size(); ++i) {
        cudaMemcpy(&d_signature_seq[i * SIGNATURE_MAX_LEN], signatures[i].seq.c_str(), signatures[i].seq.size(), cudaMemcpyHostToDevice);
    }

    // Use cudaMallocHost() to write straight to host, since only a small number of matches should be found
    cudaMallocHost(&match_scores, sizeof(double) * samples.size() * signatures.size()); // Max possible number of matches
    cudaMallocHost(&match_idx, sizeof(unsigned short) * 2 * samples.size() * signatures.size());
    // Test out zero-copy memory
    // cudaHostAlloc(&match_scores, sizeof(double) * samples.size() * signatures.size(), cudaHostAllocMapped); // Max possible number of matches
    // cudaHostAlloc(&match_idx, sizeof(unsigned short) * 2 * samples.size() * signatures.size(), cudaHostAllocMapped);
    

    // 1 Sample per block; Each thread in block corresponds to 1 signature
    getMatches<<<samples.size(), signatures.size()>>>(d_sample_seq, d_sample_qual, d_signature_seq, match_scores, match_idx, match_count);

    int h_match_count; // Number of matches after getMatches() is done
    cudaMemcpy(&h_match_count, match_count, sizeof(int), cudaMemcpyDeviceToHost);
    cudaFree(match_count);
    cudaFree(d_sample_seq);
    cudaFree(d_sample_qual);
    cudaFree(d_signature_seq);

    for (int i = 0; i < h_match_count; ++i) {
        matches.emplace_back(MatchResult(samples[match_idx[i << 1]].name,
                                         signatures[match_idx[(i << 1) + 1]].name,
                                         match_scores[i]));
    }
    
    cudaFreeHost(match_scores);
    cudaFreeHost(match_idx);
}
