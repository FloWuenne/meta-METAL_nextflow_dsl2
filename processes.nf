process run_metal {
    tag "${meta_file.baseName}"

    publishDir "${params.outdir}/${meta_file.baseName}", mode: 'copy'

    input:
    each path(meta_file)
    path sumstat_files

    output:
    val(meta_file.baseName), emit: meta_file_name
    path("METAANALYSIS1.TBL.gz"), emit: metal_out

    script:
    """
    metal ${meta_file}
    gzip METAANALYSIS1.TBL
    """
}


process plot_results {
    tag "plots.${meta_file}"
    publishDir "${params.outdir}/${meta_file}", mode: 'copy'
    
    input:
    val meta_file
    path(metal_result)

    output:
    path ("${meta_file}.meta_analysis.manhattan_plot.png")
    path ("${meta_file}.meta_analysis.qqplot_plot.png")
    path ("${meta_file}.meta_analysis.top_variants.tsv")

    script:
    """
    plot_metal_res.R --metal_res=${metal_result} --output_tag=${meta_file} --n_top_var=${params.top_n_var} --width=${params.plot_width} --height=${params.plot_height}
    """

}
/*
process extract_top_hits {
    extract_top_hits
}*/