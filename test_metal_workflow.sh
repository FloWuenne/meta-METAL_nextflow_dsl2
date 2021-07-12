## Test with different input directories for sumstat files
work_dir="$PWD"
nextflow run main.nf -resume  \
-with-report "../nextflow_reports/test_report.html" \
-with-timeline "../nextflow_reports/test_timeline.html" \
--outdir "./test_out" \
--index "./test_data/metal_list.csv" \
--plot_types "qqman" \
--vep_cache_dir "/home/florian/Postdoc/References/VEP_cache" \
--metal_scheme "SAMPLESIZE"
