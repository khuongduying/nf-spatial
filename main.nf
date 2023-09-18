#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ST } from './workflows/spatialtranscriptomics'

workflow {
    ST ()
}