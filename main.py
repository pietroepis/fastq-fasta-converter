from Bio import SeqIO
from Bio import SeqRecord
from Bio.SeqRecord import SeqRecord

# Calculate the average value of the sequence
def avg(seq):
    return round(sum(seq) / len(seq), 2)

# Calculate the percentage of an interval length in relation to the whole read length
def calc_perc(read, interval):
    return round((interval[1] - interval[0] + 1) / len(read), 2)

# Condition 1
# The read is greather than or equal to L1 and less than or equal to L2
def check_length(read, min, max):
    return len(read) >= min and len(read) <= max

# Condition 2
# The minimum quality of the read is at least equal to Q1
def check_min_quality(read, threshold):
    return min(read) > threshold

# Condition 3
# Returns the initial and final indexes of the longest subsequence in which minimum quality is greater than or equal to Q2
# and which is at least P% of the read length long
# If there isn't any interval that satisfies these two conditions, (-1, -1) is returned
def get_subregion_min_quality(read, threshold, length_perc):
    intervals = []
    start = -1

    for index in range(0, len(read)):
        if read[index] >= threshold:
            # If the previous element is less than Q2, it means that a new valid interval starts at index
            if index == 0 or read[index - 1] < threshold:
                start = index
            # If the next element is less than Q2, it means that the current open interval finishes at index
            if index == len(read) - 1 or read[index + 1] < threshold:
                intervals.append((start, index))
    
    # The list comprehension generates a list of all intervals length, the highest is taken
    max_length = max([end - start + 1 for (start, end) in intervals])

    # Only intervals that produced the max_length are taken in consideration
    max_intervals = [(start, end) for (start, end) in intervals if end - start + 1 == max_length]

    # The first interval with maximum length is returned (only if it's lenght is at least P% of thorough read length)
    return max_intervals[0] if len(max_intervals) != 0 and calc_perc(read, max_intervals[0]) >= length_perc else (-1, -1)

def write_fasta(records):
    SeqIO.write(records, "output.fa", "fasta")

L1 = int(input("L1: "))
L2 = int(input("L2: ")) 
Q1 = float(input("Q1: ")) 
Q2 = float(input("Q2: ")) 
P = float(input("P: ")) 

# L2 must be greater than L1 and Q2 must be greater than Q1
if (L2 <= L1 or Q2 <= Q1):
    print("Invalid input parameters")
    exit()

fastq_records = list(SeqIO.parse("input.fq", "fastq"))
fasta_records = []

for record in fastq_records:    
    valid_length = check_length(record.seq, L1, L2) 
    valid_min_quality = check_min_quality(record.letter_annotations["phred_quality"], Q1) 
    (subregion_start, subregion_end) = get_subregion_min_quality(record.letter_annotations["phred_quality"], Q2, P)
    
    # If current record satisfies all 3 conditions
    if valid_length and valid_min_quality and subregion_start != -1 and subregion_end != -1:
        annotated_sequence = SeqRecord(record.seq)
        annotated_sequence.id = record.id
        annotated_sequence.description = (
            "length:" + str(len(record.seq)) + " " +
            "min_quality:" + str(min(record.letter_annotations["phred_quality"])) + " " +
            "subregion_start:" + str(subregion_start) + " " +
            "subregion_end:" + str(subregion_end) + " " +
            "subregion_avg_quality:" + str(avg(record.letter_annotations["phred_quality"][subregion_start:subregion_end + 1]))
        )

        fasta_records.append(annotated_sequence)

write_fasta(fasta_records)
print("output.fa file successfully created")