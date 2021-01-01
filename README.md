# fastq-fasta-converter
FASTQ to FASTA file converter\
"Elementi di Bioinformatica" Lab Project - UniMiB

The script filters the reads of a `fastq` file provided as input (in relation to the conditions listed below) and converts them in a `fasta` file as output.

**Conditions**
1. The read length is greater than or equal to **L1** and less than or equal to **L2**
2. The minimum quality of the read is at least **Q1**
3. The read contains a subsequence whose minimum quality is greater than or equal to **Q2** and whose lenght is at least equal to **P**% of the whole read length

The script prompts the user these five parameters (L1, L2, Q1, Q2 and P) at the beginning of the execution. The program terminates if L2 is less than L1 or Q2 is less than Q1.

**Output**
The reads that overcame all the checks, are written in the output `fasta` including the following values as `description` (of course besides the `identifier`) :
- `length`: The length of the read
- `min_quality`: The minimum quality value in the read
- `subregion_start`: Start index of the subsequence with minimum quality Q2
- `subregion_end`: End index of the subsequence with minimum quality Q2
- `subregion_avg_quality`: Average Quality of the subsequence with minimum quality Q2
