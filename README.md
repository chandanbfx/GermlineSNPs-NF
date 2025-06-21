# GermlineSNPs-NF

A Nextflow pipeline for analyzing germline single nucleotide polymorphisms (SNPs) from paired-end FASTQ files. The pipeline performs quality control with FastQC, aligns reads to a reference genome using BWA, marks duplicates with Picard, calls variants with GATK HaplotypeCaller, and annotates variants with SnpEff, all within a Singularity container for reproducibility and portability.

## Description

The pipeline is optimized for low-memory systems (minimum 8GB RAM) and uses a Singularity container to ensure consistent execution across environments. The reference genome used is chromosome 1 (`chr1.fa`), and the pipeline includes support for the `GRCh38.105` SnpEff database for variant annotation.

## Requirements

- **Nextflow**: Version 23.10.1 or later
- **Singularity**: For containerized execution
- **System**: Minimum 8GB RAM (7.5GB available)
- **Input**:
  - Paired-end FASTQ files
  - Reference genome (e.g., `chr1.fa`)
  - Samples CSV file specifying sample IDs and FASTQ paths

## Setup

1. **Clone the Repository**:
   ```bash
   git clone git@github.com:chandanbfx/GermlineSNPs-NF.git
   cd GermlineSNPs-NF
   ```

2. **Download Reference Genome**:
   Download `chr1.fa` from Google Drive and place it in `data/ref/`:
   ```bash
   # Download from: https://drive.google.com/file/d/10ao_pvFXoG58Bpi2wusxXLQgB7BlPdW-/view?usp=sharing
   # After downloading, move to data/ref/chr1.fa
   ```

3. **Download Sample Data**:
   Download test FASTQ files for sample `SRR3474860`:
   ```bash
   wget -O data/reads/sample1_r1.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR347/000/SRR3474860/SRR3474860_1.fastq.gz
   wget -O data/reads/sample1_r2.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR347/000/SRR3474860/SRR3474860_2.fastq.gz
   ```

4. **Configure Samples CSV**:
   Create or update `config/samples.csv` to include the sample data:
   ```csv
   sample,fastq_1,fastq_2
   SRR3474860,/workspace/data/reads/sample1_r1.fastq.gz,/workspace/data/reads/sample1_r2.fastq.gz
   ```

5. **Download SnpEff Database**:
   Download the `GRCh38.105` database for variant annotation:
   ```bash
   singularity exec containers/pipeline.sif snpEff download -v GRCh38.105 -c config/snpEff.config
   ```

## Usage

Run the pipeline with:
```bash
./run_pipeline.sh
```

The script builds the Singularity container (if not already built) and executes the Nextflow workflow.

## Output

Outputs are stored in `results/`:
- `results/variants/`: VCF files from GATK
- `results/annotated/`: Annotated VCF files from SnpEff

## Configuration

- `config/nextflow.config`: Defines resource limits (e.g., memory, CPUs) and Singularity bind mounts.
- `config/snpEff.config`: Specifies the SnpEff database path (`/workspace/data/snpeff`).
- `config/samples.csv`: Lists sample IDs and FASTQ file paths.

## Acknowledgments

This pipeline was developed with inspiration from various open-source bioinformatics workflows. Special thanks to:
- The Nextflow community for providing a robust workflow management system.
- The SnpEff team for their comprehensive variant annotation tool.
- The GATK and BWA developers for their gold-standard tools in variant calling and alignment.
- Public data repositories (e.g., EBI SRA) for providing sample data (`SRR3474860`).

## License

MIT License (see `LICENSE` file for details).
