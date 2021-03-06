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
meta-analysis_files :   ${params.meta_file}
sumstat_file_path   :   ${params.sumstat_file_path}
index               :   ${params.index}
vep cache dir       :   ${params.vep_cache_dir}
num_forks           :   ${params.num_forks}
Significant thresh. :   ${params.sig_thresh}


Plotting parameters
===================================
"""


/* Include processes from saige_processes.nf */
include { run_metal; annotate_variants; plot_results } from './processes'


workflow {
   
   Channel
    .fromPath(params.index)
    .splitCsv(header:true)
    .map{ row-> tuple(row.metal_script, row.sumstats_dir) }
    .set { samples_ch }

   run_metal(samples_ch)

}


workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}