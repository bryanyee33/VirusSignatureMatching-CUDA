# Makefile

NVCC        := /usr/local/cuda/bin/nvcc
COMMON_SRCS := common.cc
HEADERS := kseq/kseq.h

CXXFLAGS := -std=c++20 -O3

OUTPUT_GEN  := gen_sample
OUTPUT_GEN_SIG  := gen_sig
OUTPUT_GPU  := matcher
OUTPUT_GPU_BONUS  := matcher-bonus

.DEFAULT_TARGET: all

all: $(OUTPUT_GEN) $(OUTPUT_GEN_SIG) $(OUTPUT_GPU)

$(OUTPUT_GPU): kernel_skeleton.cu $(COMMON_SRCS)
	$(NVCC) $(CXXFLAGS) -lineinfo -dlto -arch=native -o $@ $^

$(OUTPUT_GEN): gen_sample.cc 
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OUTPUT_GEN_SIG): gen_sig.cc 
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(OUTPUT_GEN) $(OUTPUT_GEN_SIG) $(OUTPUT_GPU) $(OUTPUT_GPU_BONUS)


# bonus: $(OUTPUT_GPU_BONUS)

# $(OUTPUT_GPU_BONUS): kernel_skeleton_bonus.cu $(COMMON_SRCS)
# 	$(NVCC) $(CXXFLAGS) -lineinfo -dlto -arch=native -o $@ $^