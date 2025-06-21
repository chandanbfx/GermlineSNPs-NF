process FASTQC {
  tag "QC-${sample}"

  input:
  tuple val(sample), path(reads)

  output:
  path "fastqc/${sample}_*_fastqc.{zip,html}", emit: reports

  script:
  """
  mkdir -p fastqc
  fastqc ${reads} -o fastqc
  """
}
