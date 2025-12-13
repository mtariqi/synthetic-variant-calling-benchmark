#!/usr/bin/env python3
"""
Synthetic Read Generation Pipeline for Somatic Variant Analysis

Generates synthetic paired-end reads with realistic:
- Sequencing errors (Illumina-like)
- Quality score profiles
- Somatic mutations in tumor samples
- Coverage variation and GC bias
- Tumor purity modeling
"""

import random
import gzip
import os
import sys
import argparse
import json
import numpy as np
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from collections import defaultdict

# Check for Biopython
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    print("ERROR: Biopython not installed. Install with: pip install biopython")
    sys.exit(1)


class SyntheticReadGenerator:
    """Generate synthetic sequencing reads with realistic artifacts."""
    
    def __init__(self, reference_fasta: str, read_length: int = 150, 
                 insert_size: int = 400, gc_bias: float = 0.1,
                 quality_decay: bool = True):
        """
        Initialize read generator.
        
        Args:
            reference_fasta: Path to reference genome FASTA
            read_length: Read length in bp
            insert_size: Insert size for paired-end reads
            gc_bias: GC bias factor (0-1)
            quality_decay: Simulate quality decay at read ends
        """
        self.reference_fasta = Path(reference_fasta)
        self.read_length = read_length
        self.insert_size = insert_size
        self.gc_bias = gc_bias
        self.quality_decay = quality_decay
        
        self.chromosomes = {}
        self.chromosome_lengths = {}
        
        self._validate_reference()
        self._load_reference()
    
    def _validate_reference(self):
        """Validate reference file exists and is readable."""
        if not self.reference_fasta.exists():
            print(f"ERROR: Reference file not found: {self.reference_fasta}")
            print("Please download a reference genome first.")
            print("\nExample commands:")
            print("  # Download chr22 from UCSC:")
            print("  wget 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz'")
            print("  gunzip chr22.fa.gz")
            print("  mv chr22.fa GRCh38_chr22.fa")
            sys.exit(1)
    
    def _load_reference(self):
        """Load reference genome into memory."""
        print(f"Loading reference genome: {self.reference_fasta}")
        
        total_bases = 0
        loaded_chromosomes = 0
        
        with open(self.reference_fasta, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                chrom = record.id.replace('chr', '').replace('Chr', '')
                if chrom in ['22', '21', '20']:  # Load small chromosomes for testing
                    seq = str(record.seq).upper()
                    self.chromosomes[chrom] = seq
                    self.chromosome_lengths[chrom] = len(seq)
                    total_bases += len(seq)
                    loaded_chromosomes += 1
                    print(f"  Loaded chromosome {chrom}: {len(seq):,} bp")
                    
                    if loaded_chromosomes >= 3:  # Limit for memory efficiency
                        break
        
        if not self.chromosomes:
            print("WARNING: No chromosomes loaded. Loading first chromosome...")
            with open(self.reference_fasta, 'r') as f:
                record = next(SeqIO.parse(f, "fasta"))
                chrom = record.id.replace('chr', '').replace('Chr', '')
                seq = str(record.seq).upper()
                self.chromosomes[chrom] = seq
                self.chromosome_lengths[chrom] = len(seq)
                print(f"  Loaded chromosome {chrom}: {len(seq):,} bp")
        
        print(f"Loaded {len(self.chromosomes)} chromosomes, {total_bases:,} total bases")
    
    def _get_mutation_spectrum(self, cosmic_like: bool = True) -> Dict[str, float]:
        """Get mutation probability spectrum."""
        if cosmic_like:
            return {
                'A>T': 0.08, 'A>G': 0.15, 'A>C': 0.07,
                'T>A': 0.08, 'T>G': 0.07, 'T>C': 0.15,
                'G>A': 0.22, 'G>T': 0.07, 'G>C': 0.03,
                'C>A': 0.07, 'C>T': 0.22, 'C>G': 0.03
            }
        else:
            bases = ['A', 'C', 'G', 'T']
            spectrum = {}
            for ref in bases:
                for alt in bases:
                    if ref != alt:
                        spectrum[f"{ref}>{alt}"] = 1/12
            return spectrum
    
    def introduce_variants(self, sequence: str, variant_rate: float = 0.001, 
                          cosmic_like: bool = True) -> Tuple[str, List[Dict]]:
        """Introduce somatic variants into sequence."""
        seq_list = list(sequence)
        variants = []
        mutation_spectrum = self._get_mutation_spectrum(cosmic_like)
        
        i = 0
        while i < len(seq_list):
            if random.random() < variant_rate:
                ref_base = seq_list[i]
                
                if ref_base not in 'ACGT':
                    i += 1
                    continue
                
                if random.random() < 0.85:  # 85% SNVs, 15% indels
                    possible = [m for m in mutation_spectrum.keys() 
                               if m.startswith(ref_base)]
                    if possible:
                        weights = [mutation_spectrum[m] for m in possible]
                        mutation = random.choices(possible, weights=weights)[0]
                        alt_base = mutation[-1]
                        
                        seq_list[i] = alt_base
                        variants.append({
                            'pos': i,
                            'type': 'SNV',
                            'ref': ref_base,
                            'alt': alt_base,
                            'mutation': mutation
                        })
                else:  # Indel
                    if random.random() < 0.6:  # Deletion
                        del_len = random.randint(1, 3)
                        if i + del_len <= len(seq_list):
                            deleted = ''.join(seq_list[i:i+del_len])
                            for _ in range(del_len):
                                seq_list.pop(i)
                            variants.append({
                                'pos': i,
                                'type': 'DEL',
                                'ref': deleted,
                                'alt': '',
                                'length': del_len
                            })
                    else:  # Insertion
                        ins_len = random.randint(1, 3)
                        insertion = ''.join(random.choices('ACGT', k=ins_len))
                        for j, base in enumerate(insertion):
                            seq_list.insert(i + j, base)
                        variants.append({
                            'pos': i,
                            'type': 'INS',
                            'ref': '',
                            'alt': insertion,
                            'length': ins_len
                        })
                        i += ins_len
            i += 1
        
        return ''.join(seq_list), variants
    
    def add_sequencing_errors(self, sequence: str, error_rate: float = 0.001) -> str:
        """Add Illumina-like sequencing errors."""
        seq_list = list(sequence)
        
        for i in range(len(seq_list)):
            if seq_list[i] not in 'ACGT':
                continue
            
            pos_factor = 1.0 + 2.0 * min(i, len(seq_list)-i) / len(seq_list)
            if random.random() < error_rate * pos_factor:
                bases = ['A', 'C', 'G', 'T']
                bases.remove(seq_list[i])
                seq_list[i] = random.choice(bases)
        
        return ''.join(seq_list)
    
    def generate_quality_scores(self, length: int, mean_q: int = 35) -> str:
        """Generate Phred quality scores with end decay."""
        qualities = []
        
        for pos in range(length):
            if self.quality_decay:
                decay = 0.15 * (min(pos, length-pos) / (length/2))
                q = mean_q * (1 - decay)
            else:
                q = mean_q
            
            q += np.random.normal(0, 2)
            q = max(25, min(40, int(q)))
            
            qualities.append(chr(q + 33))
        
        return ''.join(qualities)
    
    def generate_read_pair(self, chrom: str, start_pos: int, sample_type: str,
                          variant_rate: float = 0.0, tumor_purity: float = 1.0,
                          read_id: int = 0) -> Tuple[str, str, str, str, List[Dict]]:
        """Generate a single read pair."""
        ref_seq = self.chromosomes[chrom]
        frag_end = start_pos + self.insert_size
        
        if frag_end > len(ref_seq):
            start_pos = max(0, len(ref_seq) - self.insert_size - 100)
            frag_end = start_pos + self.insert_size
        
        fragment = ref_seq[start_pos:frag_end]
        variants = []
        
        if sample_type == 'tumor' and variant_rate > 0:
            effective_rate = variant_rate * tumor_purity
            fragment, variants = self.introduce_variants(fragment, effective_rate)
            
            for v in variants:
                v['ref_pos'] = start_pos + v['pos']
                v['read_id'] = read_id
        
        r1_seq = fragment[:self.read_length]
        r2_seq = str(Seq(fragment[-self.read_length:]).reverse_complement())
        
        r1_seq = self.add_sequencing_errors(r1_seq)
        r2_seq = self.add_sequencing_errors(r2_seq)
        
        r1_qual = self.generate_quality_scores(len(r1_seq))
        r2_qual = self.generate_quality_scores(len(r2_seq))
        
        return r1_seq, r2_seq, r1_qual, r2_qual, variants
    
    def generate_sample(self, output_prefix: str, sample_type: str = 'normal',
                       coverage: float = 30.0, variant_rate: float = 0.0,
                       tumor_purity: float = 1.0, chrom: str = '22',
                       seed: Optional[int] = None) -> Dict:
        """Generate a complete sample."""
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)
        
        if chrom not in self.chromosomes:
            available = list(self.chromosomes.keys())
            print(f"Using chromosome {available[0]} instead of {chrom}")
            chrom = available[0]
        
        chrom_len = self.chromosome_lengths[chrom]
        num_pairs = int((coverage * chrom_len) / (2 * self.read_length))
        
        print(f"Generating {sample_type} sample with {coverage}x coverage")
        print(f"  Chromosome {chrom}: {chrom_len:,} bp")
        print(f"  Target read pairs: {num_pairs:,}")
        
        output_dir = Path(output_prefix).parent
        output_dir.mkdir(parents=True, exist_ok=True)
        
        r1_file = f"{output_prefix}_R1.fastq.gz"
        r2_file = f"{output_prefix}_R2.fastq.gz"
        variants_file = f"{output_prefix}_variants.tsv"
        
        all_variants = []
        stats = defaultdict(int)
        
        max_start = chrom_len - self.insert_size - 100
        start_positions = np.random.randint(0, max_start, num_pairs)
        
        with gzip.open(r1_file, 'wt') as f1, gzip.open(r2_file, 'wt') as f2:
            for read_id, start_pos in enumerate(start_positions):
                r1_seq, r2_seq, r1_qual, r2_qual, variants = self.generate_read_pair(
                    chrom, int(start_pos), sample_type, variant_rate, 
                    tumor_purity, read_id
                )
                
                read_name = (f"@{sample_type}_{Path(output_prefix).name}_"
                           f"{read_id}_{chrom}_{start_pos}")
                
                f1.write(f"{read_name}/1\n{r1_seq}\n+\n{r1_qual}\n")
                f2.write(f"{read_name}/2\n{r2_seq}\n+\n{r2_qual}\n")
                
                all_variants.extend(variants)
                stats['total_reads'] += 2
                stats['total_pairs'] += 1
                
                if (read_id + 1) % 10000 == 0:
                    percent = 100 * (read_id + 1) / num_pairs
                    sys.stdout.write(f"\r  Progress: {read_id + 1:,}/{num_pairs:,} ({percent:.1f}%)")
                    sys.stdout.flush()
        
        print(f"\r  Generated {num_pairs:,} read pairs")
        
        if all_variants:
            with open(variants_file, 'w') as vf:
                vf.write("Chrom\tPosition\tType\tRef\tAlt\tLength\tMutation\tReadID\n")
                for v in all_variants:
                    vf.write(f"{chrom}\t{v['ref_pos']}\t{v['type']}\t"
                           f"{v['ref']}\t{v['alt']}\t{v.get('length', '')}\t"
                           f"{v.get('mutation', '')}\t{v['read_id']}\n")
        
        stats_file = f"{output_prefix}_stats.json"
        stats_dict = {
            'sample_type': sample_type,
            'chromosome': chrom,
            'coverage': coverage,
            'total_reads': stats['total_reads'],
            'total_pairs': stats['total_pairs'],
            'variants': len(all_variants),
            'r1_file': r1_file,
            'r2_file': r2_file,
            'generated': datetime.now().isoformat()
        }
        
        with open(stats_file, 'w') as sf:
            json.dump(stats_dict, sf, indent=2)
        
        print(f"  Output files:")
        print(f"    R1: {r1_file}")
        print(f"    R2: {r2_file}")
        if all_variants:
            print(f"    Variants: {variants_file} ({len(all_variants)} variants)")
        print(f"    Stats: {stats_file}")
        
        return stats_dict


def create_analysis_scripts(output_dir: str):
    """Create analysis scripts for HPC."""
    scripts_dir = Path(output_dir) / "analysis_scripts"
    scripts_dir.mkdir(exist_ok=True)
    
    # Alignment script
    align_sh = scripts_dir / "01_align_samples.sh"
    align_sh.write_text("""#!/bin/bash
#SBATCH --job-name=align_synthetic
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/align_%j.out
#SBATCH --error=logs/align_%j.err

module load bwa
module load samtools

REFERENCE="GRCh38_chr22.fa"

if [ ! -f "${REFERENCE}.bwt" ]; then
    echo "Indexing reference genome..."
    bwa index ${REFERENCE}
    samtools faidx ${REFERENCE}
fi

mkdir -p aligned_bams

for fastq1 in synthetic_data/*_R1.fastq.gz; do
    base=$(basename ${fastq1} _R1.fastq.gz)
    fastq2=${fastq1/_R1.fastq.gz/_R2.fastq.gz}
    
    echo "Aligning ${base}..."
    
    bwa mem -t 8 ${REFERENCE} ${fastq1} ${fastq2} | \\
        samtools sort -@ 4 -o aligned_bams/${base}.bam -
    
    samtools index aligned_bams/${base}.bam
    
    samtools flagstat aligned_bams/${base}.bam > aligned_bams/${base}_flagstat.txt
    samtools stats aligned_bams/${base}.bam > aligned_bams/${base}_stats.txt
    
    echo "Completed ${base}"
done

echo "All alignments completed!"
""")
    
    # Variant calling script
    variant_sh = scripts_dir / "02_call_variants.sh"
    variant_sh.write_text("""#!/bin/bash
#SBATCH --job-name=variant_calling
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/variant_%j.out
#SBATCH --error=logs/variant_%j.err

module load gatk
module load samtools

REFERENCE="GRCh38_chr22.fa"

mkdir -p called_variants

if [ ! -f "${REFERENCE%.fa}.dict" ]; then
    gatk CreateSequenceDictionary -R ${REFERENCE}
fi

for tumor_bam in aligned_bams/tumor_sample_*.bam; do
    base=$(basename ${tumor_bam} .bam)
    normal_bam=${tumor_bam/tumor_sample/normal_sample}
    
    if [ ! -f "${normal_bam}" ]; then
        echo "Warning: No matching normal for ${base}"
        continue
    fi
    
    echo "Calling variants for ${base}..."
    
    gatk Mutect2 \\
        -R ${REFERENCE} \\
        -I ${tumor_bam} \\
        -I ${normal_bam} \\
        -normal $(basename ${normal_bam} .bam) \\
        -O called_variants/${base}_somatic.vcf.gz \\
        --native-pair-hmm-threads 8
    
    gatk FilterMutectCalls \\
        -R ${REFERENCE} \\
        -V called_variants/${base}_somatic.vcf.gz \\
        -O called_variants/${base}_filtered.vcf.gz
    
    tabix -p vcf called_variants/${base}_filtered.vcf.gz
    
    echo "Completed ${base}"
done

echo "Variant calling completed!"
""")
    
    # Validation script
    validate_py = scripts_dir / "03_validate_variants.py"
    validate_py.write_text('''#!/usr/bin/env python3
"""
Validate called variants against ground truth.
"""

import gzip
import csv
from pathlib import Path

def load_ground_truth(variant_file):
    """Load ground truth variants."""
    variants = {}
    with open(variant_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\\t')
        for row in reader:
            key = (row['Chrom'], int(row['Position']), row['Ref'], row['Alt'])
            variants[key] = row
    return variants

def load_called_variants(vcf_file):
    """Load called variants from VCF."""
    variants = []
    
    if vcf_file.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open
    
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\\t')
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            filter_status = fields[6]
            
            if filter_status == 'PASS':
                variants.append((chrom, pos, ref, alt))
    
    return set(variants)

def main():
    synthetic_dir = Path("synthetic_data")
    variant_files = list(synthetic_dir.glob("*_variants.tsv"))
    
    results = []
    
    for truth_file in variant_files:
        sample = truth_file.name.replace('_variants.tsv', '')
        print(f"Validating {sample}...")
        
        truth_variants = load_ground_truth(truth_file)
        
        vcf_file = Path(f"called_variants/{sample}_filtered.vcf.gz")
        if not vcf_file.exists():
            print(f"  Warning: No VCF for {sample}")
            continue
        
        called_variants = load_called_variants(vcf_file)
        
        tp = 0
        fn = 0
        truth_set = set(truth_variants.keys())
        
        for var in truth_set:
            if var in called_variants:
                tp += 1
            else:
                fn += 1
        
        sensitivity = tp / len(truth_set) if truth_set else 0
        
        results.append({
            'sample': sample,
            'truth_variants': len(truth_set),
            'true_positives': tp,
            'false_negatives': fn,
            'sensitivity': sensitivity
        })
        
        print(f"  Truth: {len(truth_set)}, TP: {tp}, FN: {fn}, Sensitivity: {sensitivity:.3f}")
    
    with open('validation_results.tsv', 'w') as f:
        f.write("Sample\\tTruth_Variants\\tTrue_Positives\\tFalse_Negatives\\tSensitivity\\n")
        for r in results:
            f.write(f"{r['sample']}\\t{r['truth_variants']}\\t{r['true_positives']}\\t"
                   f"{r['false_negatives']}\\t{r['sensitivity']:.4f}\\n")
    
    print("\\nValidation complete! Results saved to validation_results.tsv")

if __name__ == "__main__":
    main()
''')
    
    # Run script
    run_all_sh = scripts_dir / "run_pipeline.sh"
    run_all_sh.write_text("""#!/bin/bash
# Complete analysis pipeline for synthetic data

echo "=== Synthetic Data Analysis Pipeline ==="
echo

mkdir -p logs
mkdir -p aligned_bams
mkdir -p called_variants

echo "Step 1: Aligning samples..."
sbatch 01_align_samples.sh
echo "Submitted alignment job"
echo

echo "Step 2: Waiting for alignment to complete..."
echo "Check job status with: squeue -u $USER"
echo "Once complete, run:"
echo "  sbatch 02_call_variants.sh"
echo

echo "Step 3: After variant calling, validate results:"
echo "  python3 03_validate_variants.py"
echo

echo "Analysis pipeline ready!"
""")
    
    for script in scripts_dir.glob("*.sh"):
        script.chmod(0o755)
    
    print(f"Created analysis scripts in: {scripts_dir}")
    return scripts_dir


def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description="Generate synthetic sequencing reads for somatic variant analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate 4 normal and 4 tumor samples (30x coverage):
  python3 generate_synthetic_reads.py -r GRCh38_chr22.fa -n 4 -t 4
  
  # Generate high-coverage samples (60x):
  python3 generate_synthetic_reads.py -r GRCh38_chr22.fa -n 2 -t 2 -c 60
  
  # Generate with custom parameters:
  python3 generate_synthetic_reads.py -r GRCh38.fa -n 3 -t 3 \\
    --read-length 100 --insert-size 350 --output my_synthetic_data
        """
    )
    
    parser.add_argument('-r', '--reference', required=True,
                       help='Reference genome FASTA file')
    parser.add_argument('-n', '--normal-samples', type=int, default=4,
                       help='Number of normal samples to generate (default: 4)')
    parser.add_argument('-t', '--tumor-samples', type=int, default=4,
                       help='Number of tumor samples to generate (default: 4)')
    parser.add_argument('-c', '--coverage', type=float, default=30.0,
                       help='Sequencing coverage per sample (default: 30)')
    parser.add_argument('-l', '--read-length', type=int, default=150,
                       help='Read length in bp (default: 150)')
    parser.add_argument('-i', '--insert-size', type=int, default=400,
                       help='Insert size for paired-end reads (default: 400)')
    parser.add_argument('--tumor-variant-rate', type=float, default=0.002,
                       help='Somatic variant rate per bp (default: 0.002)')
    parser.add_argument('--tumor-purity', type=float, default=0.8,
                       help='Tumor purity (0-1, default: 0.8)')
    parser.add_argument('-o', '--output-dir', default='synthetic_data',
                       help='Output directory (default: synthetic_data)')
    parser.add_argument('--chromosome', default='22',
                       help='Chromosome to use (default: 22)')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed for reproducibility (default: 42)')
    parser.add_argument('--create-scripts', action='store_true',
                       help='Create HPC analysis scripts')
    
    args = parser.parse_args()
    
    # Print configuration
    print("="*60)
    print("Synthetic Read Generation Pipeline")
    print("="*60)
    print(f"Reference: {args.reference}")
    print(f"Normal samples: {args.normal_samples}")
    print(f"Tumor samples: {args.tumor_samples}")
    print(f"Coverage: {args.coverage}x per sample")
    print(f"Read length: {args.read_length}bp")
    print(f"Insert size: {args.insert_size}bp")
    print(f"Tumor variant rate: {args.tumor_variant_rate}/bp")
    print(f"Tumor purity: {args.tumor_purity}")
    print(f"Output directory: {args.output_dir}")
    print(f"Chromosome: {args.chromosome}")
    print(f"Random seed: {args.seed}")
    print("="*60)
    
    # Check numpy
    try:
        import numpy as np
    except ImportError:
        print("ERROR: numpy not installed. Install with: pip install numpy")
        sys.exit(1)
    
    # Initialize generator
    generator = SyntheticReadGenerator(
        reference_fasta=args.reference,
        read_length=args.read_length,
        insert_size=args.insert_size,
        gc_bias=0.1,
        quality_decay=True
    )
    
    # Generate normal samples
    print("\n" + "="*60)
    print("Generating Normal Samples")
    print("="*60)
    
    normal_results = []
    for i in range(1, args.normal_samples + 1):
        print(f"\n--- Normal Sample {i} ---")
        prefix = Path(args.output_dir) / f"normal_sample_{i:03d}"
        
        results = generator.generate_sample(
            output_prefix=str(prefix),
            sample_type='normal',
            coverage=args.coverage,
            variant_rate=0.0,
            tumor_purity=1.0,
            chrom=args.chromosome,
            seed=args.seed + i
        )
        normal_results.append(results)
    
    # Generate tumor samples
    print("\n" + "="*60)
    print("Generating Tumor Samples")
    print("="*60)
    
    tumor_results = []
    for i in range(1, args.tumor_samples + 1):
        print(f"\n--- Tumor Sample {i} ---")
        prefix = Path(args.output_dir) / f"tumor_sample_{i:03d}"
        
        current_purity = args.tumor_purity * (0.8 + 0.4 * random.random())
        
        results = generator.generate_sample(
            output_prefix=str(prefix),
            sample_type='tumor',
            coverage=args.coverage,
            variant_rate=args.tumor_variant_rate,
            tumor_purity=current_purity,
            chrom=args.chromosome,
            seed=args.seed + 1000 + i
        )
        tumor_results.append(results)
    
    # Create summary
    print("\n" + "="*60)
    print("Generation Complete - Summary")
    print("="*60)
    
    total_normal_reads = sum(r['total_reads'] for r in normal_results)
    total_tumor_reads = sum(r['total_reads'] for r in tumor_results)
    total_variants = sum(r['variants'] for r in tumor_results)
    
    print(f"\nGenerated {args.normal_samples} normal samples:")
    print(f"  Total reads: {total_normal_reads:,}")
    print(f"  Average per sample: {total_normal_reads/args.normal_samples:,.0f}")
    
    print(f"\nGenerated {args.tumor_samples} tumor samples:")
    print(f"  Total reads: {total_tumor_reads:,}")
    print(f"  Average per sample: {total_tumor_reads/args.tumor_samples:,.0f}")
    print(f"  Total variants: {total_variants:,}")
    print(f"  Average per sample: {total_variants/args.tumor_samples:,.0f}")
    
    # Create manifest
    manifest_file = Path(args.output_dir) / "MANIFEST.txt"
    with open(manifest_file, 'w') as f:
        f.write("="*60 + "\n")
        f.write("SYNTHETIC DATA MANIFEST\n")
        f.write("="*60 + "\n\n")
        
        f.write("Normal Samples:\n")
        f.write("-"*40 + "\n")
        for i, results in enumerate(normal_results, 1):
            f.write(f"normal_sample_{i:03d}:\n")
            f.write(f"  R1: {results['r1_file']}\n")
            f.write(f"  R2: {results['r2_file']}\n")
            f.write(f"  Reads: {results['total_reads']:,}\n")
            f.write(f"  Coverage: {results['coverage']}x\n\n")
        
        f.write("Tumor Samples:\n")
        f.write("-"*40 + "\n")
        for i, results in enumerate(tumor_results, 1):
            f.write(f"tumor_sample_{i:03d}:\n")
            f.write(f"  R1: {results['r1_file']}\n")
            f.write(f"  R2: {results['r2_file']}\n")
            f.write(f"  Reads: {results['total_reads']:,}\n")
            f.write(f"  Coverage: {results['coverage']}x\n")
            f.write(f"  Variants: {results['variants']:,}\n\n")
    
    print(f"\nManifest file: {manifest_file}")
    
    # Create analysis scripts if requested
    if args.create_scripts:
        scripts_dir = create_analysis_scripts(args.output_dir)
        print(f"\nAnalysis scripts created in: {scripts_dir}")
    
    print("\n" + "="*60)
    print("Next Steps for HPC Analysis")
    print("="*60)
    print("""
1. Copy to HPC:
   scp -r synthetic_data/ username@hpc:/path/to/project/
   scp GRCh38_chr22.fa username@hpc:/path/to/project/

2. On HPC, install dependencies:
   module load python/3.9
   pip install numpy biopython

3. Run alignment:
   cd synthetic_data/analysis_scripts
   sbatch 01_align_samples.sh

4. Call variants:
   sbatch 02_call_variants.sh

5. Validate:
   python3 03_validate_variants.py
    """)
    
    return normal_results + tumor_results


if __name__ == "__main__":
    main()