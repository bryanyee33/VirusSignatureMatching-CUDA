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
#SBATCH --time=01:00:00

## Probably no need to change anything here
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --cpu_bind=core

## May want to change this depending on how much host memory you need
#SBATCH --mem-per-cpu=20G

## Just useful logfile names
#SBATCH --output=r_a2_%j.slurmlog
#SBATCH --error=r_a2_%j.slurmlog


echo "Job is running on $(hostname), started at $(date)"

cd generated/sample_seq_len
../../gen_sig 1000 3000 10000 0.1 > sig.fasta
../../gen_sample sig.fasta 2000 20 1 2 100000 100000 10 30 0.1 > samp_100k.fastq
../../gen_sample sig.fasta 2000 20 1 2 120000 120000 10 30 0.1 > samp_120k.fastq
../../gen_sample sig.fasta 2000 20 1 2 140000 140000 10 30 0.1 > samp_140k.fastq
../../gen_sample sig.fasta 2000 20 1 2 160000 160000 10 30 0.1 > samp_160k.fastq
../../gen_sample sig.fasta 2000 20 1 2 180000 180000 10 30 0.1 > samp_180k.fastq
../../gen_sample sig.fasta 2000 20 1 2 200000 200000 10 30 0.1 > samp_200k.fastq
cd ../..

cd generated/signature_seq_len
../../gen_sig 1000 3000 3000 0.1 > sig_3k.fasta
../../gen_sig 1000 4000 4000 0.1 > sig_4k.fasta
../../gen_sig 1000 5000 5000 0.1 > sig_5k.fasta
../../gen_sig 1000 6000 6000 0.1 > sig_6k.fasta
../../gen_sig 1000 7000 7000 0.1 > sig_7k.fasta
../../gen_sig 1000 8000 8000 0.1 > sig_8k.fasta
../../gen_sig 1000 9000 9000 0.1 > sig_9k.fasta
../../gen_sig 1000 10000 10000 0.1 > sig_10k.fasta
../../gen_sample sig_3k.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_3k.fastq
../../gen_sample sig_4k.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_4k.fastq
../../gen_sample sig_5k.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_5k.fastq
../../gen_sample sig_6k.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_6k.fastq
../../gen_sample sig_7k.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_7k.fastq
../../gen_sample sig_8k.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_8k.fastq
../../gen_sample sig_9k.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_9k.fastq
../../gen_sample sig_10k.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_10k.fastq
cd ../..

cd generated/sample_count
../../gen_sig 1000 3000 10000 0.1 > sig.fasta
../../gen_sample sig.fasta 2178 22 1 2 100000 200000 10 30 0.1 > samp_2.2k.fastq
../../gen_sample sig.fasta 1984 20 1 2 100000 200000 10 30 0.1 > samp_2.0k.fastq
../../gen_sample sig.fasta 1784 18 1 2 100000 200000 10 30 0.1 > samp_1.8k.fastq
../../gen_sample sig.fasta 1584 16 1 2 100000 200000 10 30 0.1 > samp_1.6k.fastq
../../gen_sample sig.fasta 1386 14 1 2 100000 200000 10 30 0.1 > samp_1.4k.fastq
../../gen_sample sig.fasta 1188 12 1 2 100000 200000 10 30 0.1 > samp_1.2k.fastq
../../gen_sample sig.fasta 990 10 1 2 100000 200000 10 30 0.1 > samp_1.0k.fastq
cd ../..

cd generated/signature_count
../../gen_sig 500 3000 10000 0.1 > sig_500.fasta
../../gen_sig 600 3000 10000 0.1 > sig_600.fasta
../../gen_sig 700 3000 10000 0.1 > sig_700.fasta
../../gen_sig 800 3000 10000 0.1 > sig_800.fasta
../../gen_sig 900 3000 10000 0.1 > sig_900.fasta
../../gen_sig 1000 3000 10000 0.1 > sig_1000.fasta
../../gen_sample sig_500.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_500.fastq
../../gen_sample sig_600.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_600.fastq
../../gen_sample sig_700.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_700.fastq
../../gen_sample sig_800.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_800.fastq
../../gen_sample sig_900.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_900.fastq
../../gen_sample sig_1000.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_1000.fastq
cd ../..

cd generated/sample_n_ratio
../../gen_sig 1000 3000 10000 0.1 > sig.fasta
../../gen_sample sig.fasta 2000 20 1 2 100000 200000 10 30 0.0 > samp_0.0.fastq
../../gen_sample sig.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_0.1.fastq
../../gen_sample sig.fasta 2000 20 1 2 100000 200000 10 30 0.2 > samp_0.2.fastq
../../gen_sample sig.fasta 2000 20 1 2 100000 200000 10 30 0.4 > samp_0.4.fastq
../../gen_sample sig.fasta 2000 20 1 2 100000 200000 10 30 0.6 > samp_0.6.fastq
../../gen_sample sig.fasta 2000 20 1 2 100000 200000 10 30 0.8 > samp_0.8.fastq
../../gen_sample sig.fasta 2000 20 1 2 100000 200000 10 30 1.0 > samp_1.0.fastq
cd ../..

cd generated/signature_n_ratio
../../gen_sig 1000 3000 10000 0.0 > sig_0.0.fasta
../../gen_sig 1000 3000 10000 0.1 > sig_0.1.fasta
../../gen_sig 1000 3000 10000 0.2 > sig_0.2.fasta
../../gen_sig 1000 3000 10000 0.4 > sig_0.4.fasta
../../gen_sig 1000 3000 10000 0.6 > sig_0.6.fasta
../../gen_sig 1000 3000 10000 0.8 > sig_0.8.fasta
../../gen_sig 1000 3000 10000 1.0 > sig_1.0.fasta
../../gen_sample sig_0.0.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_0.0.fastq
../../gen_sample sig_0.1.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_0.1.fastq
../../gen_sample sig_0.2.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_0.2.fastq
../../gen_sample sig_0.4.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_0.4.fastq
../../gen_sample sig_0.6.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_0.6.fastq
../../gen_sample sig_0.8.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_0.8.fastq
../../gen_sample sig_1.0.fasta 2000 20 1 2 100000 200000 10 30 0.1 > samp_1.0.fastq
cd ../..

cd generated/sample_virus_percentage
../../gen_sig 1000 3000 10000 0.1 > sig.fasta
../../gen_sample sig.fasta 2200 0 1 2 100000 200000 10 30 0.1 > samp_0.0.fastq
../../gen_sample sig.fasta 1980 220 1 2 100000 200000 10 30 0.1 > samp_0.1.fastq
../../gen_sample sig.fasta 1760 440 1 2 100000 200000 10 30 0.1 > samp_0.2.fastq
../../gen_sample sig.fasta 1320 880 1 2 100000 200000 10 30 0.1 > samp_0.4.fastq
../../gen_sample sig.fasta 880 1320 1 2 100000 200000 10 30 0.1 > samp_0.6.fastq
../../gen_sample sig.fasta 440 1760 1 2 100000 200000 10 30 0.1 > samp_0.8.fastq
../../gen_sample sig.fasta 0 2200 1 2 100000 200000 10 30 0.1 > samp_1.0.fastq
cd ../..

cd generated/sample_virus_count
../../gen_sig 1000 3000 10000 0.1 > sig.fasta
../../gen_sample sig.fasta 2000 20 1 1 100000 200000 10 30 0.1 > samp_1.fastq
../../gen_sample sig.fasta 2000 20 500 500 100000 200000 10 30 0.1 > samp_500.fastq
../../gen_sample sig.fasta 2000 20 1000 1000 100000 200000 10 30 0.1 > samp_1000.fastq
../../gen_sample sig.fasta 2000 20 2000 2000 100000 200000 10 30 0.1 > samp_2000.fastq
../../gen_sample sig.fasta 2000 20 3000 3000 100000 200000 10 30 0.1 > samp_3000.fastq
../../gen_sample sig.fasta 2000 20 6000 6000 100000 200000 10 30 0.1 > samp_6000.fastq
../../gen_sample sig.fasta 2000 20 8000 8000 100000 200000 10 30 0.1 > samp_8000.fastq
../../gen_sample sig.fasta 2000 20 10000 10000 100000 200000 10 30 0.1 > samp_10000.fastq
cd ../..


echo -e "\n====> Finished running.\n"
echo -e "\nJob completed at $(date)"
