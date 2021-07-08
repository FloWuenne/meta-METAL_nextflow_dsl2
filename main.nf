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
meta-analysis_files   :  ${params.meta_file}
sumstat_files         :  ${params.sumstat_files}


Plotting parameters
===================================
"""


/* Include processes from saige_processes.nf */
include { run_metal; plot_results } from './processes'


workflow {

   Channel
      .fromPath(params.meta_file)
      .ifEmpty { exit 1, "Cannot find METAL script file : ${params.meta_file}" } 
      .set{ meta_file_ch }

   Channel
      .fromPath(params.sumstat_files)
      .ifEmpty { exit 1, "Cannot find summary statistics files : ${params.sumstat_files}" }
      .collect() 
      .set{ sumstat_dir_ch }

   sumstat_dir_ch.view()

   run_metal(meta_file_ch, sumstat_dir_ch)

   plot_results(run_metal.out.meta_file_name,run_metal.out.metal_out)

    /*extract_top_hits() */


}


workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}