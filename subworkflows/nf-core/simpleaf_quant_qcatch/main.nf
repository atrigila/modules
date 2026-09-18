include { QCATCH         } from '../../../modules/nf-core/qcatch'
include { SIMPLEAF_INDEX } from '../../../modules/nf-core/simpleaf/index'
include { SIMPLEAF_QUANT } from '../../../modules/nf-core/simpleaf/quant'

workflow SIMPLEAF_QUANT_QCATCH {
    take:
    ch_genome_fasta     // (mandatory) [ meta, genome_fasta ]
    ch_genome_gtf       // (mandatory) [ meta, genome_gtf ]
    transcript_fasta    // (optional)  [ transcript_fasta ]
    simpleaf_index      // (optional)  [ simpleaf_index ]
    txp2gene            // (optional)  [ txp2gene ]
    barcode_whitelist   // (optional)  [ barcode_whitelist ]
    chemistry           // (mandatory) [ chemistry ]
    qcatch_chemistry    // (optional)  [ qcatch_chemistry ]
    skip_qcatch         // (mandatory) [ skip_qcatch ]
    resolution          // (mandatory) [ resolution ]
    ch_fastq            // (optional)  [ meta, [ fastq ] ]
    map_dir             // (optional)  [ map_dir ]

    main:
    if (!simpleaf_index && !map_dir) {
        if (transcript_fasta) {
            ch_genome_fasta_gtf = channel.of([[:], [], []])
            ch_transcript_fasta = channel.of([[id: transcript_fasta.baseName], transcript_fasta])
            if (!txp2gene) {
                error "txp2gene is required when building an index from transcript_fasta"
            }
        } else {
            ch_genome_fasta_gtf = ch_genome_fasta.combine(ch_genome_gtf).map { meta1, fasta, _meta2, gtf -> [meta1, fasta, gtf] }
            ch_transcript_fasta = channel.of([[:], []])
        }

        SIMPLEAF_INDEX(
            ch_genome_fasta_gtf,
            ch_transcript_fasta,
            channel.of([[:], []]),
            channel.of([[:], []])
        )
        ch_simpleaf_index = SIMPLEAF_INDEX.out.index.collect()
        if (txp2gene) {
            ch_txp2gene = txp2gene
        } else {
            ch_txp2gene = SIMPLEAF_INDEX.out.t2g.collect().map { _meta, t2g -> t2g }
        }
    } else if (simpleaf_index) {
        ch_simpleaf_index = simpleaf_index
        ch_txp2gene = txp2gene ?: channel.empty()
    } else {
        ch_simpleaf_index = channel.of([[:], []])
        ch_txp2gene = txp2gene ? channel.of(txp2gene) : channel.empty()
    }

    if (map_dir) {
        ch_chemistry_reads = channel.of([[:], [], []])
        ch_index_t2g = channel.of([[:], [], []])
        ch_map_dir = map_dir.map { directory -> [[id: directory.baseName], directory] }
    } else {
        ch_chemistry_reads = ch_fastq.map { meta, reads -> [meta + [chemistry: chemistry], chemistry, reads] }
        ch_index_t2g = ch_txp2gene ? ch_simpleaf_index.combine(ch_txp2gene).collect() : ch_simpleaf_index.map { meta, index -> [meta, index, []] }.collect()
        ch_map_dir = channel.of([[:], []])
    }

    SIMPLEAF_QUANT(
        ch_chemistry_reads,
        ch_index_t2g,
        channel.of([[:], barcode_whitelist ? 'unfiltered-pl' : 'knee', [], barcode_whitelist ?: []]),
        resolution,
        ch_map_dir
    )


    ch_qcatch_report = channel.empty()
    if (!skip_qcatch) {
        QCATCH(SIMPLEAF_QUANT.out.quant.map { meta, quant_dir -> [meta, qcatch_chemistry, quant_dir] })
        ch_qcatch_report = QCATCH.out.report
    }

    emit:
    txp2gene        = ch_txp2gene                                   // [ meta, txp2gene ]
    index           = ch_simpleaf_index                             // [ meta, index ]
    map             = map_dir ? ch_map_dir : SIMPLEAF_QUANT.out.map // [ meta, map ]
    quant           = SIMPLEAF_QUANT.out.quant                      // [ meta, quant ]
    qcatch_report   = ch_qcatch_report                              // [ meta, qcatch_report ]
}
