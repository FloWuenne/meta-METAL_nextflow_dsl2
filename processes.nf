process run_metal {
    tag "${meta_file.baseName}"

    publishDir "${params.outdir}/${phenoFile.baseName}/SAIGE_out_step1", mode: 'copy'

    input:
    val meta_file

    output:
    val(meta_file.baseName), emit: meta_file
    tuple path("METAANALYSIS1.TBL"), path("METAANALYSIS1.TBL.info"), emit: metal_out

    script:
    """
    metal ${params.meta_file}
    """

}

process reformat_results{

}


process plot_results {

}

process extract_top_hits {
    extract_top_hits
}