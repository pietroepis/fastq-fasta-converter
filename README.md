# fastq-fasta-converter
FASTQ to FASTA file converter\
"Elementi di Bioinformatica" Lab Project - UniMiB

The script filters the reads of a `fastq` file provided as input (in relation to the conditions listed below) and converts them in a `fasta` file as output.

**Conditions**
1. The read length is greater than or equal to **L1** and less than or equal to **L2**
2. The minimum quality of the read is at least **Q1**
3. The read contains a subsequence whose minimum quality is greater than or equal to **Q2** and whose lenght is at least equal to **P**% of the whole read length
