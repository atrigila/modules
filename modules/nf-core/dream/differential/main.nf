process DREAM_DIFFERENTIAL {
    tag "${meta.id} - ${meta.contrast_id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d6/d6fa8a7908dd357484d63bb574a378f1739440a2c0752ce91b1e0e9d1ac1638c/data' :
        'community.wave.seqera.io/library/bioconductor-edger_bioconductor-variancepartition_r-optparse:ba778938d72f30c5' }"

    input:
    tuple val(meta), path(samplesheet), path(counts)

    output:
    tuple val(meta), path("*.dream.results.tsv")        , emit: results
    tuple val(meta), path("*.dream.mean_difference.png"), emit: md_plot
    tuple val(meta), path("*.MArrayMM.dream.rds")       , emit: rdata
    tuple val(meta), path("*.dream.model.txt")          , emit: model
    tuple val(meta), path('*.dream.contrasts_plot.png') , emit: contrasts_png
    tuple val(meta), path("*.R_sessionInfo.log")        , emit: session_info
    tuple val(meta), path("*.normalised_counts.tsv")    , emit: normalised_counts, optional: true
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'dream_de.R'

    stub:
    prefix = task.ext.prefix   ?: "${meta.id}"
    """
    touch "${meta.contrast_id}.dream.results.tsv"
    touch "${meta.contrast_id}.dream.mean_difference.png"
    touch "${meta.contrast_id}.dream.contrasts_plot.png"
    touch "${meta.contrast_id}.MArrayMM.dream.rds"
    touch "${meta.contrast_id}.dream.model.txt"
    touch "${meta.contrast_id}.R_sessionInfo.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-edger: \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
        r-variancepartition: \$(Rscript -e "library(variancePartition); cat(as.character(packageVersion('variancePartition')))")
    END_VERSIONS
    """
}
