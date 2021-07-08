#!/usr/bin/env nextflow

/* Author : Florian Wuennemann */

/* Enable DSL2 syntax */
nextflow.enable.dsl=2

/* log info */
log.info """\
NEXTFLOW - DSL2 - Meta-analysis using METAL with nextflow
===================================

Meta-analysis parameters
===================================
grm_plink_input   :  ${params.grm_plink_input}


Plotting parameters
===================================
bgen_prefix       :  ${params.bgen_prefix}
bgen_suffix       :  ${params.bgen_suffix}
bgen_path         :  ${params.bgen_path}
chromosomes       :  ${params.chrom}
sampleFile        :  ${params.sampleFile}
vcfField          :  ${params.vcfField}
minMAC            :  ${params.minMAC}
minMAF            :  ${params.minMAF}
"""

/* Include processes from saige_processes.nf */
include {  } from './processes'


workflow {

   Channel
      .fromPath(params.meta_file)
      .ifEmpty { exit 1, "Cannot find pheno_File file : ${params.meta_file}" } 
      .set{ meta_file_ch }

    run_metal(meta_file_ch)

    /*reformat_results()

    plot_results()

    extract_top_hits() */


}


workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}