process GATK4_CHECKREFERENCECOMPATIBILITY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(reference_to_compare, stageAs: 'refcomp?/*', arity: '1..*'), path(reference_to_compare_fai, stageAs: 'refcomp?/*'), path(reference_to_compare_dict, stageAs: 'refcomp?/*')

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_name = input.name.toLowerCase()
    def input_flag = (input_name.endsWith('.bam') || input_name.endsWith('.cram'))
        ? '--input'
        : (input_name.endsWith('.vcf') || input_name.endsWith('.vcf.gz'))
            ? '--variant'
            : null

    if (!input_flag) {
        throw new IllegalArgumentException("Unsupported input file '${input.name}'. Expected an input ending in .bam, .cram, .vcf, or .vcf.gz.")
    }

    def reference_compare_args = reference_to_compare.collect { "--references-to-compare ${it}" }.join(' ')

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK CheckReferenceCompatibility] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        CheckReferenceCompatibility \\
        ${input_flag} ${input} \\
        ${reference_compare_args} \\
        --output ${prefix}.tsv \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
