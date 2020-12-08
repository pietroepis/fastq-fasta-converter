from Bio import SeqIO
from Bio import SeqRecord
from Bio.SeqRecord import SeqRecord

def avg(seq):
    return round(sum(seq) / len(seq), 2)

def check_length(read, min, max):
    return len(read) >= min and len(read) <= max

def check_min_quality(read, threshold):
    return min(read) > threshold

def get_subregion_min_quality(read, threshold, length_perc):
    start, end = -1, -1

    for index in range(0, len(read) - 1):
        if read[index] >= threshold:
            if index == 0 or read[index - 1] < threshold:
                start = index
            if index == len(read) - 2 or read[index + 1] < threshold:
                end = index
                if (end - start + 1) / len(read) >= length_perc:
                    break

    if (end - start + 1) / len(read) >= length_perc:
        return (start, end)
    else:
        return (-1, -1)

L1 = 30 # int(input("L1: "))
L2 = 60 # int(input("L2: "))
Q1 = 30 # float(input("Q1: "))
Q2 = 63 # float(input("Q2: "))
P = 0.15 # float(input("P: "))

fastq_records = list(SeqIO.parse("input.fq", "fastq"))
fasta_records = []

for record in fastq_records:
    valid_length = check_length(record.seq, L1, L2) 
    valid_min_quality = check_min_quality(record.letter_annotations["phred_quality"], Q1) 
    (subregion_start, subregion_end) = get_subregion_min_quality(record.letter_annotations["phred_quality"], Q2, P)

    annotated_sequence = SeqRecord(record.seq)
    annotated_sequence.id = record.id
    annotated_sequence.annotations = {
        "length": len(record.seq),
        "min_quality": Q1,
        "subregion_start":  subregion_start,
        "subregion_end":  subregion_end,
        "subregion_avg_quality": avg(record.letter_annotations["phred_quality"][subregion_start:subregion_end + 1])
    }

    if valid_length and valid_min_quality and subregion_start != -1 and subregion_end != -1 :
        fasta_records.append(annotated_sequence)