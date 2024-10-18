#include "kseq/kseq.h"
#include "common.h"

#define SAMPLE_MAX_LEN 200000
#define SIGNATURE_MAX_LEN 10000

__global__ void getMatches(const char* __restrict d_sample_seq, const char* __restrict d_sample_qual, const char* __restrict d_signature_seq,
        double* __restrict match_scores) {

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

                match_scores[blockIdx.x * blockDim.x + threadIdx.x] = (double)tot / j;
                return;

            } else if (!sample[i + j] || sample[i + j] != 'N' && signature[j] != 'N' &&
                    sample[i + j] != signature[j]) { // No match
                break;
            }
        }
    }
}

void runMatcher(const std::vector<klibpp::KSeq>& samples, const std::vector<klibpp::KSeq>& signatures, std::vector<MatchResult>& matches) {
    cudaStream_t stream1, stream2, stream3;
    cudaStreamCreate(&stream1);
    cudaStreamCreate(&stream2);
    cudaStreamCreate(&stream3);

    int sample_arr_len = samples.size() * SAMPLE_MAX_LEN;
    int signature_arr_len = signatures.size() * SIGNATURE_MAX_LEN;

    // Device variables
    char *d_sample_seq, *d_sample_qual, *d_signature_seq;
    int *match_count; // 1 int storing the number of matches
    cudaMallocAsync(&d_sample_seq, sample_arr_len, stream1);
    cudaMemsetAsync(d_sample_seq, 0, sample_arr_len, stream1); // Initialise arrays to '\0'

    cudaMallocAsync(&d_sample_qual, sample_arr_len, stream2);
    cudaMemsetAsync(d_sample_qual, 0, sample_arr_len, stream2);

    cudaMallocAsync(&d_signature_seq, signature_arr_len, stream3);
    cudaMemsetAsync(d_signature_seq, 0, signature_arr_len, stream3);

    cudaMallocAsync(&match_count, sizeof(int), stream3);
    cudaMemsetAsync(match_count, 0, 1, stream3); // Initialise match_count to 0

    // Pinned memory
    double *match_scores;      // [Score_1, Score2, ...]

    for (int i = 0; i < signatures.size(); ++i) {
        cudaMemcpyAsync(&d_signature_seq[i * SIGNATURE_MAX_LEN], signatures[i].seq.c_str(), signatures[i].seq.size(), cudaMemcpyHostToDevice, stream3);
    }
    for (int i = 0; i < samples.size(); ++i) {
        cudaMemcpyAsync(&d_sample_seq[i * SAMPLE_MAX_LEN], samples[i].seq.c_str(), samples[i].seq.size(), cudaMemcpyHostToDevice, stream1);
        cudaMemcpyAsync(&d_sample_qual[i * SAMPLE_MAX_LEN], samples[i].qual.c_str(), samples[i].qual.size(), cudaMemcpyHostToDevice, stream2);
    }
    

    // Use cudaMallocHost() to write straight to host, since only a small number of matches should be found
    cudaMallocHost(&match_scores, sizeof(double) * samples.size() * signatures.size()); // Max possible number of matches
    memset(match_scores, 0, sizeof(double) * samples.size() * signatures.size());

    // 1 Sample per block; Each thread in block corresponds to 1 signature
    getMatches<<<samples.size(), signatures.size()>>>(d_sample_seq, d_sample_qual, d_signature_seq, match_scores);

    cudaFreeAsync(d_sample_seq, stream1);
    cudaFreeAsync(d_sample_qual, stream2);
    cudaFreeAsync(d_signature_seq, stream3);
    cudaDeviceSynchronize();
    for (int i = 0; i < samples.size() * signatures.size(); ++i) {
        if (match_scores[i] != 0) {
            matches.emplace_back(MatchResult(samples[i / signatures.size()].name,
                                             signatures[i % signatures.size()].name,
                                             match_scores[i]));
        }
    }
    
    cudaFreeHost(match_scores);
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
    cudaStreamDestroy(stream3);
}
