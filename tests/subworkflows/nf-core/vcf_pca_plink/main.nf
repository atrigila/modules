#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_PCA_PLINK } from '../../../../subworkflows/nf-core/vcf_pca_plink/main.nf'

workflow test_vcf_pca_plink {


    input = Channel.of([
        [id:'sampleID1', chr:'22'],
        file("https://github.com/nf-core/test-datasets/raw/gwas/data/data_shrink_chunk_4500/chr22.vcf.bgz", checkIfExists: true),
        file("https://github.com/nf-core/test-datasets/raw/gwas/data/data_shrink_chunk_4500/chr22.vcf.bgz.tbi", checkIfExists: true)
    ])

    fasta = file("https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta", checkIfExists: true)
    fasta_fai = file("https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai", checkIfExists: true)

    VCF_PCA_PLINK ( input, fasta )
}
