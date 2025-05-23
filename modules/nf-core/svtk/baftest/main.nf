process SVTK_BAFTEST {
    tag "$meta.id"
    label 'process_single'

    //Conda is not supported at the moment: https://github.com/broadinstitute/gatk-sv/issues/787
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svtk:0.0.20190615--py37h73a75cf_2':
        'biocontainers/svtk:0.0.20190615--py37h73a75cf_2' }"

    input:
    tuple val(meta), path(bed), path(baf), path(baf_index), path(batch)

    output:
    tuple val(meta), path("*.metrics")  , emit: metrics
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def VERSION = '0.0.20190615' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    svtk baf-test \\
        ${bed} \\
        ${baf} \\
        --batch ${batch} \\
        ${args} \\
        > ${prefix}.metrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svtk: ${VERSION}
    END_VERSIONS
    """
}
