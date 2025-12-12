#!/usr/bin/env python3
"""
Generate synthetic paired-end FASTQ reads for:
  - NORMAL sample
  - TUMOR sample (with somatic mutations)

Creates:
  normal_R1.fastq.gz
  normal_R2.fastq.gz
  tumor_R1.fastq.gz
  tumor_R2.fastq.gz

Total reads for each = 5000

Author: Md Tariqul Islam (Tariq)
"""

import gzip
import random
from Bio import SeqIO

# -----------------------------
# PARAMETERS
# -----------------------------
reference_fasta = "reference.fasta"
chromosome = "chr1"
start = 100000
end   = 100600   # 600 bp region gives room for fragments

num_reads = 5000           # Number of paired reads per sample
read_length = 150
mutation_rate_tumor = 0.02  # Tumor only — 2% mutations
insert_mean = 350
insert_std  = 30

random_seed = 42
random.seed(random_seed)

bases = ["A", "C", "G", "T"]

# -----------------------------
# MUTATION FUNCTIONS
# -----------------------------
def introduce_snp(base):
    return random.choice([b for b in bases if b != base])

def mutate_sequence(seq, rate):
    seq_list = list(seq)
    positions = []
    for i, b in enumerate(seq_list):
        if random.random() < rate:
            seq_list[i] = introduce_snp(b)
            positions.append(i)
    return "".join(seq_list), positions

def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def random_qual(n):
    return "".join(chr(random.randint(33, 73)) for _ in range(n))


# -----------------------------
# LOAD REFERENCE REGION
# -----------------------------
print("[INFO] Loading reference...")
ref = SeqIO.to_dict(SeqIO.parse(reference_fasta, "fasta"))

if chromosome not in ref:
    raise ValueError(f"Chromosome {chromosome} not found in reference.")

ref_seq = str(ref[chromosome].seq[start:end])
print(f"[INFO] Region length: {len(ref_seq)} bp")


# -----------------------------
# CREATE NORMAL SAMPLE (NO MUTATIONS)
# -----------------------------
print("\n[INFO] Generating NORMAL reads...")

normal_R1 = gzip.open("normal_R1.fastq.gz", "wt")
normal_R2 = gzip.open("normal_R2.fastq.gz", "wt")

for i in range(1, num_reads + 1):

    insert_len = max(read_length * 2, int(random.gauss(insert_mean, insert_std)))
    pos = random.randint(0, len(ref_seq) - insert_len)

    fragment = ref_seq[pos:pos + insert_len]
    seq_r1 = fragment[:read_length]
    seq_r2 = reverse_complement(fragment[-read_length:])

    q1 = random_qual(read_length)
    q2 = random_qual(read_length)

    read_id = f"@NORMAL_read_{i}"

    normal_R1.write(f"{read_id}/1\n{seq_r1}\n+\n{q1}\n")
    normal_R2.write(f"{read_id}/2\n{seq_r2}\n+\n{q2}\n")

normal_R1.close()
normal_R2.close()

print("[INFO] NORMAL sample generated.")


# -----------------------------
# CREATE TUMOR SAMPLE (WITH SOMATIC MUTATIONS)
# -----------------------------
print("\n[INFO] Generating TUMOR reads...")

# Mutate the reference for tumor
tumor_ref, tumor_mut_positions = mutate_sequence(ref_seq, mutation_rate_tumor)
print(f"[INFO] Tumor mutations introduced: {len(tumor_mut_positions)}")
print("Mutation positions (example):", tumor_mut_positions[:10])

tumor_R1 = gzip.open("tumor_R1.fastq.gz", "wt")
tumor_R2 = gzip.open("tumor_R2.fastq.gz", "wt")

for i in range(1, num_reads + 1):

    insert_len = max(read_length * 2, int(random.gauss(insert_mean, insert_std)))
    pos = random.randint(0, len(tumor_ref) - insert_len)

    fragment = tumor_ref[pos:pos + insert_len]
    seq_r1 = fragment[:read_length]
    seq_r2 = reverse_complement(fragment[-read_length:])

    q1 = random_qual(read_length)
    q2 = random_qual(read_length)

    read_id = f"@TUMOR_read_{i}"

    tumor_R1.write(f"{read_id}/1\n{seq_r1}\n+\n{q1}\n")
    tumor_R2.write(f"{read_id}/2\n{seq_r2}\n+\n{q2}\n")

tumor_R1.close()
tumor_R2.close()

print("[INFO] TUMOR sample generated.")


print("\n[INFO] Synthetic dataset creation complete:")
print("  normal_R1.fastq.gz")
print("  normal_R2.fastq.gz")
print("  tumor_R1.fastq.gz")
print("  tumor_R2.fastq.gz")
print("\nReady for BWA → BAM → Variant Calling → hap.py benchmarking.")

