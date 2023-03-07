// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { BCFTOOLS_NORM   as BCFTOOLS_NORM_VCF   } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM   as BCFTOOLS_NORM_BCF   } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_ANNOTATE                      } from '../../../modules/nf-core/bcftools/annotate/main'

include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'

workflow VCF_PCA_PLINK {

    take:
    // TODO nf-core: edit input (take) channels
    ch_vcf          // channel: [ val(meta), vcf, tbi ]
    ch_reference    // channel: [ val(meta), fasta , fai ]

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    BCFTOOLS_NORM_VCF ( ch_vcf )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM_VCF.out.versions.first())

    BCFTOOLS_ANNOTATE ( BCFTOOLS_NORM.out.vcf )
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())

    BCFTOOLS_NORM_BCF( BCFTOOLS_ANNOTATE.out.vcf )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM_BCF.out.versions.first())

    BCFTOOLS_INDEX( BCFTOOLS_ANNOTATE.out.bcf )
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())


    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

