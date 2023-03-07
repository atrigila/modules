#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_PCA_PLINK } from '../../../../subworkflows/nf-core/vcf_pca_plink/main.nf'

workflow test_vcf_pca_plink {
    
    input = file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)

    VCF_PCA_PLINK ( input )
}
