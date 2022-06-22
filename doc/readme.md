# Potato Mapping Tests

## Summary

The files included here are in pairs:

1) The general TSV stats
2) the Pathogen counts TSV

## Tests

### stats_Potato_blast

1) Trim with porechop --middle_threshold 75 --require_two_barcodes
2) create combined_db from both St_path_db.fasta and Potato sequences
3) Mapped to combined_db using blastn -task megablast -word_size 12 -evalue 1e-80 -outfmt 6
4) Stats gathered using: cut -d $'\t' -f 1 | sort | uniq | wc -l


### stats_Potato_bwa

1) Trim with porechop --middle_threshold 75 --require_two_barcodes
2) Mapped trimmed fastq file to Potato sequences using BWA -x ont2d
3) Extract fastq unmapped sequences
4) Mapped to St_path_db.fasta using BWA -x ont2d
5) Stats gathered from BAM files using: cut -d $'\t' -f 1 | sort | uniq | wc -l

### stats_Potato_bwa-blast-no_word

1) Trim with porechop --middle_threshold 75 --require_two_barcodes
2) Mapped trimmed fastq file to Potato sequences using BWA -x ont2d
3) Extract fastq unmapped sequences
4) Mapped to St_path_db.fasta using blastn -task megablast -evalue 1e-80 -perc_identity 90 -outfmt 6
5) Stats gathered using: cut -d $'\t' -f 1 | sort | uniq | wc -l

### stats_Potato_bwa-blast-word_size

1) Trim with porechop --middle_threshold 75 --require_two_barcodes
2) Mapped trimmed fastq file to Potato sequences using BWA -x ont2d
3) Extract fastq unmapped sequences
4) Mapped to St_path_db.fasta using blastn -task megablast -word_size 12 -evalue 1e-80 -perc_identity 90 -outfmt 6
5) Stats gathered using: cut -d $'\t' -f 1 | sort | uniq | wc -l

### stats_Potato_bwa-graphmap-85

1) Trim with porechop --middle_threshold 75 --require_two_barcodes
2) Mapped trimmed fastq file to Potato sequences using BWA -x ont2d
3) Extract fastq unmapped sequences
4) Mapped to St_path_db.fasta using graphmap align -z 1e-85
5) Stats gathered from BAM files using: cut -d $'\t' -f 1 | sort | uniq | wc -l

- Stats are also filtered using a python script that filters by percent identity 90
- filter_identity.py uses the MD:Z Field in BAM files with equation:
  - identity = 100 * M / (M+S)
  - Where M are matches and S are SNPs

### stats_Potato_bwa-graphmap-90

1) Trim with porechop --middle_threshold 75 --require_two_barcodes
2) Mapped trimmed fastq file to Potato sequences using BWA -x ont2d
3) Extract fastq unmapped sequences
4) Mapped to St_path_db.fasta using graphmap align -z 1e-90
5) Stats gathered from BAM files using: cut -d $'\t' -f 1 | sort | uniq | wc -l

- Stats are also filtered using a python script that filters by percent identity 90
- filter_identity.py uses the MD:Z Field in BAM files with equation:
  - identity = 100 * M / (M+S)
  - Where M are matches and S are SNPs

### stats_Potato_bwa-graphmap-95

1) Trim with porechop --middle_threshold 75 --require_two_barcodes
2) Mapped trimmed fastq file to Potato sequences using BWA -x ont2d
3) Extract fastq unmapped sequences
4) Mapped to St_path_db.fasta using graphmap align -z 1e-95
5) Stats gathered from BAM files using: cut -d $'\t' -f 1 | sort | uniq | wc -l

- Stats are also filtered using a python script that filters by percent identity 90
- filter_identity.py uses the MD:Z Field in BAM files with equation:
  - identity = 100 * M / (M+S)
  - Where M are matches and S are SNPs
