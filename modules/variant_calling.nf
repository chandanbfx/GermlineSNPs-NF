process CREATE_SEQUENCE_DICTIONARY {
    input:
    path ref

    output:
    path "${ref.baseName}.dict"

    script:
    """
    gatk CreateSequenceDictionary -R ${ref} -O ${ref.baseName}.dict
    """
}

process GATK_HAPLOTYPECALLER {
    tag "VC-${sample}"
    publishDir "${params.outdir}/vcf", mode: 'copy' // Save output to results/vcf

    input:
    path bam
    path ref
    path index_files // Expects chr1.fa.fai and chr1.dict
    val sample

    output:
    path "${bam.baseName}.g.vcf.gz", emit: vcf

    script:
    """
    gatk --java-options "-Xmx6G" HaplotypeCaller \
        -I ${bam} \
        -R ${ref} \
        -O ${bam.baseName}.g.vcf.gz \
        --sample-name ${sample}
    """
}
