#!/bin/bash
#PBS -o 98_log_files/Qsub_trinity.out
#PBS -l walltime=180:00:00
#PBS -l mem=160g
#PBS -l ncpus=12
#PBS -q omp

###PBS -r n

# Move to present working directory
#cd $PBS_O_WORKDIR

#TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
#SCRIPT=$0
#NAME=$(basename $0)
#LOG_FOLDER="98_log_files"
#cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"


# Load module
. /appli/bioinfo/trinity/2.8.5/env.sh


# Global variables
DATA_DIR="/home1/datawork/tdestanq/dev_larve_2019/Transcriptome_Saccostrea/results/trimming/"			#path to directory containing trimmed fastq files
WORKING_DIRECTORY="/home1/scratch/tdestanq/Qsub_Trinity/"

# Variable
SEQ_TYPE="fq"
LEFT="${DATA_DIR}/NEBNext_dual_i5_B9.T3_1_trim_R1_paired.fastq.gz,${DATA_DIR}/NEBNext_dual_i5_D8.T1_1_trim_R1_paired.fastq.gz,${DATA_DIR}/NEBNext_dual_i5_E9.T4_1_trim_R1_paired.fastq.gz,${DATA_DIR}/NEBNext_dual_i5_G8.T2_1_trim_R1_paired.fastq.gz"
RIGHT="${DATA_DIR}/NEBNext_dual_i5_B9.T3_1_trim_R2_paired.fastq.gz,${DATA_DIR}/NEBNext_dual_i5_D8.T1_1_trim_R2_paired.fastq.gz,${DATA_DIR}/NEBNext_dual_i5_E9.T4_1_trim_R2_paired.fastq.gz,${DATA_DIR}/NEBNext_dual_i5_G8.T2_1_trim_R2_paired.fastq.gz"
NB_CPU=12
MAX_MEM="150G"
ADD_PARAMS=""			#Specify additionnal parameters for Trinity if needed

#Creating working directory if not existing
mkdir -p $WORKING_DIRECTORY


#Running Trinity

time Trinity --seqType fq --max_memory 150 --left $LEFT --right $RIGHT --CPU 12 $ADD_PARAMS --output $WORKING_DIRECTORY >& 98_log_files/Qsub_Trinity.log 2>&1


