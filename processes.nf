process run_metal {
    tag "${meta_file.baseName}"

    publishDir "${params.outdir}/${meta_file.baseName}", mode: 'copy'

    input:
    tuple path(meta_file), path(sumstat_file_path)
    

    output:
    val(meta_file.baseName), emit: meta_file_name
    path("METAANALYSIS1.TBL"), emit: metal_out
    path("METAANALYSIS1.TBL.info"), emit: metal_out_info

    script:
    """
    ls ${sumstat_file_path} > all_files.txt
    for file_name in `cut -d "/" -f2 all_files.txt`
    do 
        ln -s ${sumstat_file_path}/\${file_name} ./\${file_name}
    done
    metal ${meta_file}
    """
}

process annotate_variants {
    tag "${meta_file}"
    publishDir "${params.outdir}/${meta_file}", mode: 'copy'

    input:
    val meta_file
    path(metal_result)
    path(vep_cache_dir)
    val(metal_scheme)

    output:
    path ("sumstats_reformat.tsv")
    path ("${meta_file}.FUMA_format.full.tsv.gz")
    path ("${meta_file}.FUMA_format.final.tsv.gz")
    path ("${meta_file}.VEP_annotation.tsv")
    path ("${meta_file}.VEP_annotation.tsv_summary.html")
    path ("${meta_file}.sig_HET.tsv.gz")

    script:
     if( metal_scheme == 'SAMPLESIZE')
        """
        ## Sort Meta-analysis results for VEP
        mkdir tmp_dir
        awk 'gsub(/(:| )+/,"\t")' ${metal_result} | sort -t" " -Vk1,2 -T ./tmp_dir > sorted_metal_out.tsv

        ## VEP
        awk -F " " '{print \$1"\t"\$2"\t"\$2"\t"\$3"/"\$4"\t+\t"\$1"_"\$2"_"\$3"_"\$4}' sorted_metal_out.tsv > sumstats_reformat.tsv
        vep -i sumstats_reformat.tsv -o ${meta_file}.VEP_annotation.tsv --offline --species homo_sapiens --dir_cache ${vep_cache_dir}  --tab --pick --check_existing --symbol --fork ${params.num_forks}

        ## Format for FUMA online annotation : https://fuma.ctglab.nl
        echo "rsID\tID\tSYMBOL\tIMPACT" > rsIDs.txt 
        grep -v "#" ${meta_file}.VEP_annotation.tsv | awk -F " " '{print \$13"\t"\$1"\t"\$18"\t"\$14}' >> rsIDs.txt 

        echo "CHR\tPOS\tREF\tALT\tPVAL\tZSCORE" > sumstats_metal_sub.tsv
        awk -F " " '{print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$13"\t"\$12}' sorted_metal_out.tsv >> sumstats_metal_sub.tsv

        ## Merge rsID information with meta-analysis parameters for FUMA formatted output
        paste -d "\t" rsIDs.txt sumstats_metal_sub.tsv > ${meta_file}.FUMA_format.full.tsv
        awk -F " " '\$1 != "-" {print \$0}' ${meta_file}.FUMA_format.full.tsv >> ${meta_file}.FUMA_format.final.tsv
        gzip ${meta_file}.FUMA_format.full.tsv
        gzip ${meta_file}.FUMA_format.final.tsv 

        ## Get all variants with significant heterogeneity
        head -1 sorted_metal_out.tsv > ${meta_file}.sig_HET.tsv
        tail -n +2 sorted_metal_out.tsv | awk -F " " '\$15 < 0.05 {print \$0}' >> ${meta_file}.sig_HET.tsv
        gzip ${meta_file}.sig_HET.tsv
        """

    else if(metal_scheme == 'STDERR')
        """
        ## Sort Meta-analysis results for VEP
        mkdir tmp_dir
        awk 'gsub(/(:| )+/,"\t")' ${metal_result} | sort -t" " -Vk1,2 -T ./tmp_dir > sorted_metal_out.tsv

        ## VEP
        awk -F " " '{print \$1"\t"\$2"\t"\$2"\t"\$3"/"\$4"\t+\t"\$1"_"\$2"_"\$3"_"\$4}' sorted_metal_out.tsv > sumstats_reformat.tsv
        vep -i sumstats_reformat.tsv -o ${meta_file}.VEP_annotation.tsv --offline --species homo_sapiens --dir_cache ${vep_cache_dir}  --tab --pick --check_existing --symbol --fork ${params.num_forks}

        ## Format for FUMA online annotation: https://fuma.ctglab.nl
        echo "rsID\tID\tSYMBOL\tIMPACT" > rsIDs.txt
        grep -v "#" ${meta_file}.VEP_annotation.tsv | awk -F " " '{print \$13"\t"\$1"\t"\$18"\t"\$14}' >> rsIDs.txt

        echo "CHR\tPOS\tREF\tALT\tPVAL\tBETA\tSE"  > sumstats_metal_sub.tsv
        awk -F " " '{print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$13"\t"\$11"\t"\$12}' sorted_metal_out.tsv  >> sumstats_metal_sub.tsv

        ## Merge rsID information with meta-analysis parameters for FUMA formatted output
        paste -d "\t" rsIDs.txt sumstats_metal_sub.tsv > ${meta_file}.FUMA_format.full.tsv
        awk -F " " '\$1 != "-" {print \$0}' ${meta_file}.FUMA_format.full.tsv >> ${meta_file}.FUMA_format.final.tsv
        gzip ${meta_file}.FUMA_format.full.tsv
        gzip ${meta_file}.FUMA_format.final.tsv

        ## Get all variants with significant heterogeneity
        head -1 sorted_metal_out.tsv > ${meta_file}.sig_HET.tsv
        tail -n +2 sorted_metal_out.tsv | awk -F " " '\$15 < 0.05 {print \$0}' >> ${meta_file}.sig_HET.tsv
        gzip ${meta_file}.sig_HET.tsv
        """

}


process plot_results {
    tag "${meta_file}"
    publishDir "${params.outdir}/${meta_file}", mode: 'copy'
    
    input:
    val meta_file
    path(metal_result)

    output:
    path ("${meta_file}.meta_analysis.manhattan_plot.png")
    path ("${meta_file}.meta_analysis.qqplot_plot.png")
    path ("${meta_file}.meta_analysis.sig_variants.tsv")

    script:
    """
    plot_metal_res.R --metal_res=${metal_result} --output_tag=${meta_file} --sig_thresh=${params.sig_thresh} --width=${params.plot_width} --height=${params.plot_height} --plot_types=${params.plot_types}
    """
}