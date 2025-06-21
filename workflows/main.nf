println "DEBUG: params.ref_genome = ${params.ref_genome}"
println "DEBUG: params.samples = ${params.samples}"
if (!file(params.samples).exists()) {
    error "ERROR: Samples file not found: ${params.samples}"
}
if (!params.ref_genome) {
    error "ERROR: params.ref_genome is not defined"
}
ref = file(params.ref_genome)
ref_fai = Channel.fromPath("${params.ref_genome}.fai")
params.outdir = "results"

Channel.fromPath(params.samples)
  | splitCsv(header: true)
  | map { row -> [ row.sample, [ file(row.fastq_1), file(row.fastq_2) ] ] }
  | set { read_pairs }

include { FASTQC } from '../modules/qc'
include { MARKDUPLICATES } from '../modules/markdups'
include { GENERATE_BWA_INDEX } from '../modules/alignment'
include { ALIGNMENT } from '../modules/alignment'
include { CREATE_SEQUENCE_DICTIONARY; GATK_HAPLOTYPECALLER } from '../modules/variant_calling'
include { ANNOTATE } from '../modules/annotation'

params.test = false
if (params.test) {
    params.reads = "data/test/*_{1,2}.fastq.gz"
    params.ref_genome = "data/test/chr1.fa"
    params.snpeff_db = "data/test/snpeff_db"
    println "RUNNING IN TEST MODE"
}

workflow {
    GENERATE_BWA_INDEX(ref)
    bwa_indexes = GENERATE_BWA_INDEX.out.collect()
    CREATE_SEQUENCE_DICTIONARY(ref)
    ref_indexes = ref_fai.mix(CREATE_SEQUENCE_DICTIONARY.out).collect()
    FASTQC(read_pairs)
    ALIGNMENT(read_pairs, ref, bwa_indexes)
    MARKDUPLICATES(ALIGNMENT.out)
    // Combine sample name with BAM file
    bam_with_sample = read_pairs.combine(MARKDUPLICATES.out.bam)
        | map { sample, fastqs, bam -> [sample, bam] }
    bam_with_sample.view { "DEBUG: bam_with_sample = $it" }
    GATK_HAPLOTYPECALLER(
        bam_with_sample.map { it[1] }, // bam
        ref,
        ref_indexes,
        bam_with_sample.map { it[0] } // sample
    )
    ANNOTATE(GATK_HAPLOTYPECALLER.out, "GRCh38.105")
}
