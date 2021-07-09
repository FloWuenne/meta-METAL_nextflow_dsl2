process run_metal {
    tag "${meta_file.baseName}"

    publishDir "${params.outdir}/${meta_file.baseName}", mode: 'copy'

    input:
    tuple path(meta_file), path(sumstat_files)
    

    output:
    val(meta_file.baseName), emit: meta_file_name
    path("${sumstat_files}/METAANALYSIS1.TBL"), emit: metal_out

    script:
    """
    cp ${meta_file} ${sumstat_files}
    cd ${sumstat_files}
    metal ${meta_file}
    """
}

process annotate_variants {
    tag "plots.${meta_file}"
    publishDir "${params.outdir}/${meta_file}", mode: 'copy'

    input:
    val meta_file
    path(metal_result)

    output:
    path ("${meta_file}.FUMA_format.tsv.gz")
    path ("${meta_file}.VEP_annotation.tsv.gz")


    script:
    """
    ## VEP
    awk 'gsub(/(:| )+/,"\t")' ${metal_result} | gawk -F " " '{print $1"\t"$2"\t"$2"\t"$3"/"$4"\t+\t"$1"_"$2"_"$3"_"$4}' > sumstats_reformat.tsv
    vep -i sumstats_reformat.tsv -o ${meta_file}.VEP_annotation.tsv --offline --cache --dir_cache ${params.vep_cache_dir} --species homo_sapiens --tab --pick

    ## Format for FUMA online annotation : https://fuma.ctglab.nl/
    echo "rsID/tPVAL/tCHR\tPOS\tREF\tALT\tBETA\tSE"  > sumstats.FUMA_format.tsv
    awk 'gsub(/(:| )+/,"\t")' ${metal_result} | gawk -F " " '{print }' >> sumstats.FUMA_format.tsv
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
    plot_metal_res.R --metal_res=${metal_result} --output_tag=${meta_file} --n_top_var=${params.top_n_var} --width=${params.plot_width} --height=${params.plot_height} --plot_types=${params.plot_types}
    """

}