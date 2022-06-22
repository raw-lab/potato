#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=Potato_BLASTn
#SBATCH --nodes=4
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100GB
#SBATCH --time=10-0
#SBATCH -o slurm-%x-%j.out
#SBATCH --mail-type=END,FAIL,REQUEUE

########## README ##########
# potato_sbatch.sh
# This script utilizes SLURM on a cluster to speed up processing.

# porechop   0.2.4
# blast      2.11.0+
# seqkit     0.16.1
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
module load blast
module load seqkit
# Porechop is in an Anaconda environment since it was not pre-installed in my cluster.
module load anaconda3
eval "$(conda shell.bash hook)"
conda activate /projects/raw_lab/envs/aligners


# CPUs per node to use with BWA/Graphmap
CPUS=$SLURM_CPUS_ON_NODE

# Location of the database files
DB_BLAST_POTATO_PATHOGEN=database/combined_db

echo Indexing Databases
makeblastdb -in database/combined.fasta -dbtype nucl -title Combined -out $DB_BLAST_POTATO_PATHOGEN


EVAL="1e-80"
IDENTITY="90"

# main output directory
OUTPATH=$SLURM_JOB_NAME
mkdir -p $OUTPATH

# location of the stats file to write to
STATS=$OUTPATH/stats_$SLURM_JOB_NAME.txt

# List of files to process
FILE_LIST=(samples/*.fastq)

echo "====================================================="
echo "Trimming Files                              : $(date +%H:%M)"
echo "====================================================="
for FILE in ${FILE_LIST[@]}
do
    # basename of the sample file
    name=$(basename $FILE .fastq)
    # results folder for the sample
    RESULTS=$OUTPATH/$name
    mkdir -p $RESULTS
    # output trimmed file
    TRIMMED=$RESULTS/$name-trimmed.fastq

    srun --nodes=1 --ntasks=1 \
        porechop -t $CPUS --middle_threshold 75 --require_two_barcodes -i $FILE -o $TRIMMED > $RESULTS/poreshop.log &
done
wait # wait for all jobs to finish before moving on


echo "====================================================="
echo "BLASTn to pathogen DB using e-value: $EVAL : $(date +'%D %T')"
echo "====================================================="
for FILE in ${FILE_LIST[@]}
do
    name=$(basename $FILE .fastq)
    RESULTS=$OUTPATH/$name
    TRIMMED=$RESULTS/$name-trimmed.fastq
    # output tsv
    BLAST_OUT=$RESULTS/pathogen_blast.tsv    

    # start blastn on a node,
    # save the TSV file (for stats later),
    echo -n > $BLAST_OUT
    srun --nodes=1 --ntasks=1 --quiet \
        cat $TRIMMED | seqkit fq2fa | \
        blastn -num_threads $CPUS -task megablast -word_size 12 -evalue $EVAL -outfmt 6 -db $DB_BLAST_POTATO_PATHOGEN >> $BLAST_OUT &
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
    TRIMMED=$RESULTS/$name-trimmed.fastq
    BLAST_OUT=$RESULTS/pathogen_blast.tsv

    count=$(awk '{s++}END{print s/4}' $TRIMMED)
    potato_reads=0
    pathogen_match=$(cat $BLAST_OUT | cut -d $'\t' -f 1 | sort | uniq | wc -l)
    unmapped_reads=$(( count - pathogen_match ))

    unmatched=$(( count - pathogen_match ))
    echo Sequences in $name: $count >> $STATS
    echo Mapped potato reads: $potato_reads >> $STATS
    echo Unmapped potato reads: $unmapped_reads >> $STATS
    echo Mapped pathogen reads: $pathogen_match >> $STATS
    echo Unmapped reads: $unmatched >> $STATS

    echo "Matched pathogen counts: $name" >> $STATS
    while read line; do
        abv=${line:0:7}
        count=$(cat $BLAST_OUT | grep $abv | cut -d $'\t' -f 1 | sort | uniq | wc -l)
        echo "  $abv:$count" >> $STATS
    done <pathogen_list.txt

done

./parse_stats.py $STATS $OUTPATH/$(basename $STATS .txt).tsv

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "Ran in $SECONDS seconds"
echo "======================================================"
echo ""

