process MARKDUPLICATES {
  tag "DEDUP-${sample}"

  input:
  path bam

  output:
  path "${bam.baseName}.dedup.bam", emit: bam
  path "${bam.baseName}.dedup.metrics", emit: metrics

  script:
  """
     gatk MarkDuplicates \
    -I ${bam} \
    -O ${bam.baseName}.dedup.bam \
    -M ${bam.baseName}.dedup.metrics
  samtools index ${bam.baseName}.dedup.bam
  """
}
