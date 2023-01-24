#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SOMALIER_EXTRACT  } from '../../../../../modules/nf-core/somalier/extract/main.nf'
include { SOMALIER_ANCESTRY } from '../../../../../modules/nf-core/somalier/ancestry/main.nf'
include { UNTAR             } from '../../../../../modules/nf-core/untar/main.nf'

workflow test_somalier_ancestry {


    ancestry_labels_1kg = file("https://github.com/brentp/somalier/raw/73db124d3fe9febe3a53787707554f863595b48f/scripts/ancestry-labels-1kg.tsv", checkIfExists: true)
    1kg_sommalier = file("https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz", checkIfExists: true)
    whole_genome_sommalier = file("https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/somalier/whole_genome.somalier", checkIfExists: true)
    
    labels  = [[id:"1kg"], ancestry_labels_1kg]
    labelled_somalier_tar = [[id:"1kg"], 1kg_sommalier]
    input = [[id:"test"], whole_genome_sommalier]

    UNTAR ( labelled_somalier_tar )
    ch_labelled_somalier_files = labels.join(UNTAR.out.untar)

    SOMALIER_ANCESTRY ( input, ch_labelled_somalier_files )

}
