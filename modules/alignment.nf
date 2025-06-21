// Process definitions
process GENERATE_BWA_INDEX {
  
  input: path ref
  output: path("${ref}.*")  // Emits .amb, .ann, etc.

  script:
  """
  bwa index ${ref}
  """
}

process ALIGNMENT {
    
    input:
    tuple val(sample), path(reads)
    path ref
    path index_files

    output:
    path "${sample}.bam"

    script:
    """
    bwa mem ${ref} ${reads} | samtools sort -o ${sample}.bam
    """
}

