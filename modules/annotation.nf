process ANNOTATE {
  tag "ANN-${sample}"
  publishDir "results/annotated", mode: 'copy'

  input:
  path vcf
  val genome_version  // Pre-downloaded database (e.g., GRCh38.105)

  output:
  path "${vcf.baseName}.annotated.vcf", emit: vcf

  script:
  """
      snpEff -Xmx3G -noStats -v -c /workspace/config/snpEff.config ${genome_version} ${vcf} > ${vcf.baseName}.annotated.vcf
  """
}
