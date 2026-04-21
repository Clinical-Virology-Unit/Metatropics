#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/metatropics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/metatropics

    Website: https://nf-co.re/metatropics
    Slack  : https://nfcore.slack.com/channels/metatropics
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Resolve friendly host keywords (e.g. --Host human) into local cached FASTA paths
HostReferences.apply(params, log, baseDir)

def providedHumanHostFasta = params.Human_host_fasta
def providedOtherHostFasta = params.Other_host_fasta
def legacyHumanParamSupplied = params.fasta ? true : false
def legacyOtherParamSupplied = params.host_fasta ? true : false

def genomeHumanHost   = WorkflowMain.getGenomeAttribute(params, 'Human_host_fasta')
def genomeLegacyFasta = WorkflowMain.getGenomeAttribute(params, 'fasta')
def genomeOtherHost   = WorkflowMain.getGenomeAttribute(params, 'Other_host_fasta')

if (!params.Human_host_fasta) {
    params.Human_host_fasta = genomeHumanHost ?: genomeLegacyFasta ?: params.fasta
}

if (!params.Other_host_fasta) {
    params.Other_host_fasta = genomeOtherHost ?: params.host_fasta
}

if (legacyHumanParamSupplied && !providedHumanHostFasta) {
    log.warn "Parameter '--fasta' is deprecated. Please use '--Human_host_fasta' instead."
}

if (legacyOtherParamSupplied && !providedOtherHostFasta) {
    log.warn "Parameter '--host_fasta' is deprecated. Please use '--Other_host_fasta' instead."
}

// Keep legacy parameters in sync for backwards compatibility
params.fasta      = params.Human_host_fasta
params.host_fasta = params.Other_host_fasta

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { METATROPICS } from './workflows/metatropics'

//
// WORKFLOW: Run main nf-core/metatropics analysis pipeline
//
workflow NFCORE_METATROPICS {
    METATROPICS ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_METATROPICS ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
