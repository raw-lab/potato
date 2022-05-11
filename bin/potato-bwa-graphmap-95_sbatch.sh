#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=Potato_bwa-graphmap-95
#SBATCH --nodes=5
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=100GB
#SBATCH --time=1-0
#SBATCH -o slurm-%x-%j.out
#SBATCH --mail-type=END,FAIL,REQUEUE

########## README ##########
# potato_sbatch.sh
# This script utilizes SLURM on a cluster to speed up processing.

# porechop   0.2.4
# bwa        0.7.17-r1188
# graphmap   0.5.2
# samtools   1.10 (using htslib 1.10)
# python
# pandas
############################


echo "====================================================="
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""

# Track the time it takes to run the script
SECONDS=0

# Load dependencies
module load bwa
module load samtools
# Porechop and Graphmap are in an Anaconda environment since they were not pre-installed in my cluster.
module load anaconda3
eval "$(conda shell.bash hook)"
conda activate /projects/raw_lab/envs/aligners

# CPUs per node to use with BWA/Graphmap
CPUS=$SLURM_CPUS_ON_NODE

# Location of the database files
DB_POTATO=database/potato_db.fasta
DB_PATHOGEN=database/St_path_db.fasta

# evalue cutoff for mapping to pathogens
EVAL="1e-80"
IDENTITY="95"

# main output directory
OUTPATH=$SLURM_JOB_NAME
mkdir -p $OUTPATH

# location of the stats file to write to
STATS=$OUTPATH/stats_$SLURM_JOB_NAME.txt

# List of files to process
FILE_LIST=(samples/*.fastq)

echo "====================================================="
echo "Trimming Files (porechop)          : $(date +'%D %T')"
echo "====================================================="
for FILE in ${FILE_LIST[@]}
do
    # basename of the sample file
    name=$(basename $FILE .fastq)
    # results folder for the sample
    RESULTS=$OUTPATH/$name
    mkdir -p $RESULTS
    # output trimmed file
    TRIMMED=$RESULTS/potato_trimmed.fastq

    # start porechop on an available node
    srun --nodes=1 --ntasks=1 --quiet \
        porechop -t $CPUS --middle_threshold 75 --require_two_barcodes -i $FILE -o $TRIMMED > $RESULTS/poreshop.log &
done
wait # wait for all jobs to finish before moving on


echo "====================================================="
echo "BWA to Potato (Filtering)      : $(date +'%D %T')"
echo "====================================================="
for FILE in ${FILE_LIST[@]}
do
    name=$(basename $FILE .fastq)
    RESULTS=$OUTPATH/$name
    # input trimmed file
    TRIMMED=$RESULTS/potato_trimmed.fastq
    # output unmapped fastq file
    UNMAPPED=$RESULTS/potato_unmapped.fastq

    # start BWA on a node,
    # save the log file,
    # save the BAM file (for stats later),
    # and export the unmapped regions to a FASTQ file
    srun --nodes=1 --ntasks=1 --quiet \
        bwa mem -t $CPUS -x ont2d $DB_POTATO $TRIMMED 2> $RESULTS/potato.log | \
        samtools view -b - | tee $RESULTS/potato.bam | \
        samtools fastq -f 4 - > $UNMAPPED &
done
wait


echo "====================================================="
echo "Graphmap to pathogen DB using e-value: $EVAL : $(date +'%D %T')"
echo "====================================================="
for FILE in ${FILE_LIST[@]}
do
    name=$(basename $FILE .fastq)
    RESULTS=$OUTPATH/$name
    # input unmapped file
    UNMAPPED=$RESULTS/potato_unmapped.fastq
    # output bam file
    PATHOGEN=$RESULTS/pathogen.bam

    # start Graphmap on a node,
    # save the log file,
    # save the BAM file (for stats later),
    srun --nodes=1 --ntasks=1 --quiet \
        graphmap align -t $CPUS -z $EVAL -r $DB_PATHOGEN -d $UNMAPPED -o /dev/stdout 2> $RESULTS/pathogen.log | \
        samtools view -b - > $PATHOGEN &
done
wait


echo "====================================================="
echo "Gathering stats                             : $(date)"
echo "====================================================="

printf "" > $STATS
for FILE in ${FILE_LIST[@]}
do
    name=$(basename $FILE .fastq)
    RESULTS=$OUTPATH/$name
    TRIMMED=$RESULTS/potato_trimmed.fastq
    UNMAPPED=$RESULTS/potato_unmapped.fastq
    PATHOGEN=$RESULTS/pathogen.bam

    count=$(awk '{s++}END{print s/4}' $TRIMMED)
    potato_reads=$(samtools view -F 4 $RESULTS/potato.bam | cut -d $'\t' -f 1 | sort | uniq | wc -l)
    unmapped_reads=$(awk '{s++}END{print s/4}' $UNMAPPED)
    pathogen_match=$(samtools view -F 4 $PATHOGEN | ./filter_identity.py $IDENTITY | cut -d $'\t' -f 1 | sort | uniq | wc -l )
    unmatched=$(samtools view -f 4 $PATHOGEN | cut -d $'\t' -f 1 | sort | uniq | wc -l )

    echo Sequences in $name: $count >> $STATS
    echo Mapped potato reads: $potato_reads >> $STATS
    echo Unmapped potato reads: $unmapped_reads >> $STATS
    echo Mapped pathogen reads: $pathogen_match >> $STATS
    echo Unmapped reads: $unmatched >> $STATS

    echo "Matched pathogen counts: $name" >> $STATS
    total=0
    while read line; do
        abv=${line:0:7}
        if [ $pathogen_match -gt $total ]
        then
            count=$(samtools view -F 4 $PATHOGEN | ./filter_identity.py $IDENTITY | grep $abv | cut -d $'\t' -f 1 | sort | uniq | wc -l)
        else
            count=0
        fi
        echo "  $abv:$count" >> $STATS
        ((total=total+count))
    done <pathogen_list.txt

done

# Create TSV Table
./parse_stats.py $STATS $OUTPATH/$(basename $STATS .txt).tsv

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "Ran in $SECONDS seconds"
echo "======================================================"
