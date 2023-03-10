// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { BCFTOOLS_NORM   as BCFTOOLS_NORM_VCF   } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM   as BCFTOOLS_NORM_BCF   } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_ANNOTATE                      } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_INDEX                         } from '../../../modules/nf-core/bcftools/index/main'
include { PLINK_BCF                              } from '../../../modules/nf-core/plink/bcf/main'
include { PLINK_INDEP                            } from '../../../modules/nf-core/plink/indep/main'
include { PLINK_EXTRACT                          } from '../../../modules/nf-core/plink/extract/main'


workflow VCF_PCA_PLINK {

    take:
    // TODO nf-core: edit input (take) channels
    ch_vcf          // channel: [ [id:'samplename', chr:'22'], vcf, tbi ] Meta map requires 'id' and 'chr' values.
    ch_reference    // channel: [ fasta , fai ]

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

    PLINK_BCF( BCFTOOLS_ANNOTATE.out.bcf )
    ch_versions = ch_versions.mix(PLINK_BCF.out.versions.first())

    ch_bed = PLINK_BCF.out.bed
    ch_bim = PLINK_BCF.out.bim
    ch_fam = PLINK_BCF.out.fam
    ch_bed_bim_fam = ch_bed.join(ch_bim).join(ch_fam)

    ch_window_size = Channel.value(50)
    ch_variant_count = Channel.value(5)
    ch_variance_inflation_factor = Channel.value(1.5)

    PLINK_INDEP( ch_bed_bim_fam, ch_window_size, ch_variant_count, ch_variance_inflation_factor )
    ch_versions = ch_versions.mix(PLINK_INDEP.out.versions.first())

    ch_prunein = PLINK_INDEP.out.prunein
    ch_bed_bim_bam_prunein = ch_bed_bim_fam.join(ch_prunein)

    PLINK_EXTRACT( ch_bed_bim_bam_prunein )
    ch_versions = ch_versions.mix(PLINK_EXTRACT.out.versions.first())



    emit:
    // TODO nf-core: edit emitted channels
    bed      = PLINK_EXTRACT.out.bed           // channel: [ val(meta), [ bed ] ]
    bim      = PLINK_EXTRACT.out.bim          // channel: [ val(meta), [ bim ] ]
    fam      = PLINK_EXTRACT.out.fam          // channel: [ val(meta), [ fam ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

