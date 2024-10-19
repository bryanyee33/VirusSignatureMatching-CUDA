#include "kseq/kseq.h"
#include "common.h"

#define SAMPLE_MAX_LEN 200000
#define SIGNATURE_MAX_LEN 10000

// 32 bits / 3 bits = 10 char stored per int
// 200000 / 10 char = 20000 ints required
#define INT_SAMPLE_MAX_LEN 20000
#define INT_SIGNATURE_MAX_LEN 1000

__global__ void getMatches(const unsigned int* __restrict d_sample_seq, const char* __restrict d_sample_qual, const unsigned int* __restrict d_signature_seq,
        double* __restrict match_scores, unsigned short* __restrict match_idx, int* __restrict match_count) {

    const unsigned int *sample = &d_sample_seq[INT_SAMPLE_MAX_LEN * blockIdx.x];
    const unsigned int *signature = &d_signature_seq[INT_SIGNATURE_MAX_LEN * threadIdx.x];

    // **NOTE: Does not work properly since the sample characters may not be aligned to 32-bits
    for (int i = 0; i < INT_SAMPLE_MAX_LEN; ++i) {
        for (int j = 0; j < INT_SIGNATURE_MAX_LEN; ++j) {
            unsigned int curr_xor = sample[i + j] ^ signature[j];

            // Bitmask every 3rd bit i.e. 0b00100100100...100 == 0x24924924
            // Every 1 in such position indicates that there was an 'N' character in either the sample or signature
            // (positions with 'N' in both sample and signature are cancelled out by XOR)
            unsigned int temp = curr_xor & 0x24924924;

            // 0xC0000000 are the 2 leftmost bits. Mask it out for easier calculation (checked separately).
            // Repeat every encounter of the 'N' bit 2 more times to ignore the following 2 bits as well
            unsigned int result = curr_xor & ~(0xC0000000 | temp | temp >> 1 | temp >> 2);

            if (result == 0 && signature[j] & 0x40000000) { // The int matches, and signature ended. Found match.
                const char *qual = &d_sample_qual[(INT_SAMPLE_MAX_LEN * blockIdx.x + i) * 10];

                // Find how many chars in last int
                int char_count = 1;                    // Can't have 0 chars
                unsigned int temp = signature[j] >> 3; // Can't have 0 chars
                for (int k = 0; k < 9; ++k) {
                    if (temp & 7 == 7) { // Found last char indicator
                        break;
                    }
                    temp >>= 3;
                    ++char_count;
                }

                int tot = 0;
                for (int s = 0; s < j * 10 + char_count; ++s) {
                    tot += qual[s] - 33;
                }

                int idx = atomicAdd(match_count, 1); // Each idx will be unique (no need synchronisation)
                match_scores[idx] = (double)tot / (j * 10 + char_count);
                match_idx[idx << 1] = blockIdx.x; // Sample idx
                match_idx[(idx << 1) + 1] = threadIdx.x; // Signature idx
                return;

            } else if (result != 0 || sample[i + j] & 0x40000000) { // No match, or sample ended.
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

    
    int sample_arr_len = samples.size() * INT_SAMPLE_MAX_LEN;
    int signature_arr_len = signatures.size() * INT_SIGNATURE_MAX_LEN;

    // calloc() to initialise all to '\0'; set to max possible length for each type.
    char *h_sample_qual = (char*)calloc(samples.size() * SAMPLE_MAX_LEN, sizeof(char));
    unsigned int *h_sample_seq = (unsigned int*)malloc(sample_arr_len * sizeof(unsigned int));
    unsigned int *h_signature_seq = (unsigned int*)malloc(signature_arr_len * sizeof(unsigned int));

    // Device variables
    char *d_sample_qual;
    unsigned int *d_sample_seq, *d_signature_seq;
    int *match_count; // 1 int storing the number of matches
    cudaMallocAsync(&d_sample_qual, samples.size() * SAMPLE_MAX_LEN, stream2);
    cudaMallocAsync(&d_sample_seq, sample_arr_len, stream1);
    cudaMallocAsync(&d_signature_seq, signature_arr_len, stream3);
    cudaMallocAsync(&match_count, sizeof(int), stream3);
    cudaMemsetAsync(match_count, 0, 1, stream3); // Initialise match_count to 0

    // Pinned memory
    double *match_scores;      // [Score_1, Score2, ...]
    unsigned short *match_idx; // [Samp_idx_1, Sig_idx_1, Samp_idx_2, Sig_idx_2, ...]

    // Combine each string array into 1 for transfer
    for (int i = 0; i < samples.size(); ++i) {
        memcpy(&h_sample_qual[i * SAMPLE_MAX_LEN], &samples[i].qual[0], samples[i].qual.size());
    }
    cudaMemcpyAsync(d_sample_qual, h_sample_qual, sample_arr_len, cudaMemcpyHostToDevice, stream2);

    for (int i = 0; i < signatures.size(); ++i) {
        // Number of ints needed for this signature -1 (last int after loop)
        int int_count_minus_1 = (signatures[i].seq.size() + 9) / 10 - 1;

        for (int j = 0; j < int_count_minus_1; ++j) {
            unsigned int val = 0; // Store 10 A/C/G/T/N  into this int, each taking 3 bits
            for (int k = 0; k < 10; ++k) {
                val <<= 3;
                val += ((signatures[i].seq[j * 10 + k] - 65) >> 1) % 7;
            }
            h_signature_seq[i * INT_SIGNATURE_MAX_LEN + j] = val;
        }
        // Last int for this signature
        unsigned int val = 7; // 0b111 to indicate position of last char
        for (int k = 0; k < signatures[i].seq.size() - int_count_minus_1 * 10; ++k) {
            val <<= 3;
            val += ((signatures[i].seq[int_count_minus_1 * 10 + k] - 65) >> 1) % 7;
        }
        h_signature_seq[i * INT_SIGNATURE_MAX_LEN + int_count_minus_1] = val | 0x40000000; // Set 2nd leftmost bit as 1 to indicate last int
    }
    cudaMemcpyAsync(d_signature_seq, h_signature_seq, signature_arr_len, cudaMemcpyHostToDevice, stream3);

    for (int i = 0; i < samples.size(); ++i) {
        // Number of ints needed for this sample -1 (last int after loop)
        int int_count_minus_1 = (samples[i].seq.size() + 9) / 10 - 1;

        for (int j = 0; j < int_count_minus_1; ++j) {
            unsigned int val = 0; // Store 10 A/C/G/T/N  into this int, each taking 3 bits
            for (int k = 0; k < 10; ++k) {
                val <<= 3;
                val += ((samples[i].seq[j * 10 + k] - 65) >> 1) % 7;
            }
            h_sample_seq[i * INT_SAMPLE_MAX_LEN + j] = val;
        }
        // Last int for this sample
        unsigned int val = 7; // 0b111 to indicate position of last char
        for (int k = 0; k < samples[i].seq.size() - int_count_minus_1 * 10; ++k) {
            val <<= 3;
            val += ((samples[i].seq[int_count_minus_1 * 10 + k] - 65) >> 1) % 7;
        }
        h_sample_seq[i * INT_SAMPLE_MAX_LEN + int_count_minus_1] = val | 0x40000000; // Set 2nd leftmost bit as 1 to indicate last int
    }
    cudaMemcpyAsync(d_sample_seq, h_sample_seq, sample_arr_len, cudaMemcpyHostToDevice, stream1);
    

    // Use cudaMallocHost() to write straight to host, since only a small number of matches should be found
    cudaMallocHost(&match_scores, sizeof(double) * samples.size() * signatures.size()); // Max possible number of matches
    cudaMallocHost(&match_idx, sizeof(unsigned short) * 2 * samples.size() * signatures.size());
    // Test out zero-copy memory
    // cudaHostAlloc(&match_scores, sizeof(double) * samples.size() * signatures.size(), cudaHostAllocMapped); // Max possible number of matches
    // cudaHostAlloc(&match_idx, sizeof(unsigned short) * 2 * samples.size() * signatures.size(), cudaHostAllocMapped);


    // 1 Sample per block; Each thread in block corresponds to 1 signature
    getMatches<<</*samples.size()*/ 1, signatures.size()>>>(d_sample_seq, d_sample_qual, d_signature_seq, match_scores, match_idx, match_count);
    // free() right after async getMatches() launch, since they are only needed for copy
    free(h_sample_qual);
    free(h_sample_seq);
    free(h_signature_seq);

    int h_match_count; // Number of matches after getMatches() is done
    cudaMemcpy(&h_match_count, match_count, sizeof(int), cudaMemcpyDeviceToHost);
    cudaFreeAsync(match_count, stream3);
    cudaFreeAsync(d_sample_qual, stream2);
    cudaFreeAsync(d_sample_seq, stream1);
    cudaFreeAsync(d_signature_seq, stream3);

    for (int i = 0; i < h_match_count; ++i) {
        matches.emplace_back(MatchResult(samples[match_idx[i << 1]].name,
                                         signatures[match_idx[(i << 1) + 1]].name,
                                         match_scores[i]));
    }
    
    cudaFreeHost(match_scores);
    cudaFreeHost(match_idx);
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
    cudaStreamDestroy(stream3);
}
