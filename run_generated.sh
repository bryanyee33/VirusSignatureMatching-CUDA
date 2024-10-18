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
#SBATCH --time=02:00:00

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

echo -e "\n====> Compiling...\n"
# make clean && make

cd generated

for folder in sample_virus_count #sample_seq_len signature_seq_len sample_count signature_count sample_n_ratio signature_n_ratio sample_virus_percentage sample_virus_count
do  
    echo "Folder:$folder"
    cd $folder
    {
    for sample in $(ls | grep samp)
    do
        echo "Running:$sample"
        if [ -e "sig.fasta" ]; then
            ../../matcher $sample sig.fasta
        else
            ../../matcher $sample $(echo "$sample" | sed 's/samp/sig/g' | sed 's/fastq/fasta/g')
        fi
        echo
    done
    } #> ../"$folder.result"
    echo
    cd ..
done

#cat r_a2_12345.slurmlog | grep -P "(FOR|Running|Folder)" | cut -d":" -f2 | cat

echo -e "\n====> Finished running.\n"
echo -e "\nJob completed at $(date)"
