#!/bin/bash

## Change this to a job name you want
#SBATCH --job-name=assignment2

## Change based on length of job and `sinfo` partitions available
##SBATCH --partition=gpu

## Request for a specific type of node
## Commented out for now, change if you need one
## gpu:1 ==> any gpu. For e.g., --gres=gpu:a100-40:1 gets you one of the A100 GPU shared instances

#SBATCH --constraint=xgph
#SBATCH --gpus=a100-40
##SBATCH --gpus=h100-96
##SBATCH --constraint=xgpi

## Must change this based on how long job will take. We are just expecting 30 seconds for now
#SBATCH --time=00:30:00

## Probably no need to change anything here
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --cpu_bind=core

## May want to change this depending on how much host memory you need
#SBATCH --mem-per-cpu=20G

## Just useful logfile names
#SBATCH --output=a2_%j.slurmlog
#SBATCH --error=a2_%j.slurmlog


echo "Job is running on $(hostname), started at $(date)"

# Set the nvidia compiler directory
# NVCC=/usr/local/cuda/bin/nvcc

# Check that it exists and print some version info
# [[ -f $NVCC ]] || { echo "ERROR: NVCC Compiler not found at $NVCC, exiting..."; exit 1; }
# echo "NVCC info: $($NVCC --version)"

# Actually compile the code
echo -e "\n====> Compiling...\n"
# $NVCC -arch native -O3 --std=c++17 -o sum sum.cu
make clean && make

echo -e "\n====> Running...\n"

# ./matcher-bonus samp.fastq sig.fasta
# ./matcher samp.fastq sig.fasta

# nsys profile --stats=true ./matcher samp_100%.fastq sig.fasta

# ./matcher samp_100%.fastq sig.fasta
# ./bench-a100 samp_100%.fastq sig.fasta
# ./bench-h100 samp_100%.fastq sig.fasta

# ./matcher samp_200%.fastq sig_100%.fasta
# ./bench-a100 samp_200%.fastq sig_100%.fasta
# ./bench-h100 samp_200%.fastq sig_100%.fasta

./matcher samp.fastq sig.fasta
./bench-a100 samp.fastq sig.fasta
# ./bench-h100 samp.fastq sig.fasta


# COMBINE ARRAYS
# for (( i=0; i<5; i=i+1 )) do
#     ./match_arraycombine samp.fastq sig.fasta
# done
# for (( i=0; i<5; i=i+1 )) do
#     ./match_noarraycombine samp.fastq sig.fasta
# done
# for (( i=0; i<5; i=i+1 )) do
#     ./match_arraycombine_noasync samp.fastq sig.fasta
# done
# for (( i=0; i<5; i=i+1 )) do
#     ./match_noarraycombine_noasync samp.fastq sig.fasta
# done

# for (( i=0; i<5; i=i+1 )) do
#     ./match_arraycombine samp_fullseq.fastq sig_fullseq.fasta
# done
# for (( i=0; i<5; i=i+1 )) do
#     ./match_noarraycombine samp_fullseq.fastq sig_fullseq.fasta
# done
# for (( i=0; i<5; i=i+1 )) do
#     ./match_arraycombine_noasync samp_fullseq.fastq sig_fullseq.fasta
# done
# for (( i=0; i<5; i=i+1 )) do
#     ./match_noarraycombine_noasync samp_fullseq.fastq sig_fullseq.fasta
# done


# RETRIEVE RESULT
# for (( i=0; i<5; i=i+1 )) do
#     ./match_simplest samp.fastq sig.fasta
# done
# for (( i=0; i<5; i=i+1 )) do
#     ./match_nocalc samp.fastq sig.fasta
# done
# for (( i=0; i<5; i=i+1 )) do
#     ./match_calc samp.fastq sig.fasta
# done

# for (( i=0; i<5; i=i+1 )) do
#     ./match_simplest samp_fullmatch.fastq sig_fullseq.fasta
# done
# for (( i=0; i<5; i=i+1 )) do
#     ./match_nocalc samp_fullmatch.fastq sig_fullseq.fasta
# done
# for (( i=0; i<5; i=i+1 )) do
#     ./match_calc samp_fullmatch.fastq sig_fullseq.fasta
# done


echo -e "\n====> Finished running.\n"

echo -e "\nJob completed at $(date)"
