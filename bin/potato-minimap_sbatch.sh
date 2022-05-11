#!/bin/bash

#SBATCH --partition=Andromeda
#SBATCH --job-name=Potato_Minimap
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128GB
#SBATCH --time=1-0
#SBATCH -o slurm-%x-%j.out
#SBATCH --mail-type=END,FAIL,REQUEUE

########## README ##########
# potato_sbatch.sh
# This script utilizes SLURM on a cluster to speed up processing.
# This was one of the first runs, I had not yet optimized it for multiple nodes

# porechop   0.2.4
# minimap2   2.18-r1035-dirty
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

module load samtools
module load minimap2
module load anaconda3
eval "$(conda shell.bash hook)"
conda activate /projects/raw_lab/envs/aligners


CPUS=$SLURM_CPUS_ON_NODE

POTATO_DB=database/potato_db.mmi
PATHOGEN_DB=database/St_path_db.mmi

echo "Creating index"
minimap2 -d $POTATO_DB database/potato_db.fasta
minimap2 -d $PATHOGEN_DB database/St_path_db.fasta


STATS=stats_minimap.txt
echo "Minimap Pathogen Stats" > $STATS
ONEHIT=stats_minimap_onehit.txt
echo "Minimap single-hit Stats" > $ONEHIT

FILE_LIST=(samples/*.fastq)

for f in ${FILE_LIST[@]}
do
    name=$(basename $f .fastq)
    outfolder=results_minimap/$name
    FILE=$outfolder/$name-trimmed.fastq.gz
    echo Processing $name $f
    mkdir -p $outfolder
    unmapped=$outfolder/potato_unmapped.fastq

    echo "Removing adapters"
    porechop -t $CPUS --middle_threshold 75 --require_two_barcodes -i $f -o $FILE > $outfolder/porechop.log

    echo "Mapping to potato DB"
    minimap2 -ax map-ont $POTATO_DB $FILE -t $CPUS | \
        samtools view -b - | \
        tee $outfolder/potato.bam | \
        samtools fastq -f 4 - > $unmapped


    echo "Mapping to pathogen DB"
    minimap2 -ax map-ont $PATHOGEN_DB $unmapped -t $CPUS | samtools view -b - > $outfolder/pathogen.bam

    echo "Gathering stats"
    count=$(awk '{s++}END{print s/4}' $f)
    potato_reads=$(samtools view -F 4 $outfolder/potato.bam | cut -d $'\t' -f 1 | sort | uniq | wc -l)
    unmapped_reads=$(awk '{s++}END{print s/4}' $unmapped)
    pathogen_match=$(samtools view -F 4 $outfolder/pathogen.bam | cut -d $'\t' -f 1 | sort | uniq | wc -l )
    unmatched=$(samtools view -f 4 $outfolder/pathogen.bam | cut -d $'\t' -f 1 | sort | uniq | wc -l )
    echo Sequences in $name: $count >> $STATS
    echo Mapped potato reads: $potato_reads >> $STATS
    echo Unmapped potato reads: $unmapped_reads >> $STATS
    echo Mapped pathogen reads: $pathogen_match >> $STATS
    echo Unmapped reads: $unmatched >> $STATS

    echo "Matched pathogen counts:" >> $STATS
    while read line; do
        abv=${line:0:7}
        count=$(samtools view -F 4 $outfolder/pathogen.bam | grep $abv | cut -d $'\t' -f 1 | sort | uniq | wc -l)
        echo "  $abv:$count" >> $STATS
    done <pathogen_list.txt


    echo "Gathering One-hit Match stats"
    count=$(awk '{s++}END{print s/4}' $f)
    potato_reads=$(samtools view -F 4 $outfolder/potato.bam | cut -d $'\t' -f 1 | sort | uniq -u | wc -l)
    unmapped_reads=$(awk '{s++}END{print s/4}' $unmapped)
    pathogen_match=$(samtools view -F 4 $outfolder/pathogen.bam | cut -d $'\t' -f 1 | sort | uniq -u | wc -l )
    unmatched=$(samtools view -f 4 $outfolder/pathogen.bam | cut -d $'\t' -f 1 | sort | uniq -u | wc -l )
    echo Sequences in $name: $count >> $ONEHIT
    echo Mapped potato reads: $potato_reads >> $ONEHIT
    echo Unmapped potato reads: $unmapped_reads >> $ONEHIT
    echo Mapped pathogen reads: $pathogen_match >> $ONEHIT
    echo Unmapped reads: $unmatched >> $ONEHIT

    echo "Matched pathogen counts:" >> $ONEHIT
    while read line; do
        abv=${line:0:7}
        count=$(samtools view -F 4 $outfolder/pathogen.bam | grep $abv | cut -d $'\t' -f 1 | sort | uniq -u | wc -l)
        echo "  $abv:$count" >> $ONEHIT
    done <pathogen_list.txt

done

python parse_stats.py $STATS $(echo $STATS | cut -f 1 -d '.').tsv
python parse_stats.py $ONEHIT $(echo $ONEHIT | cut -f 1 -d '.').tsv


echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
