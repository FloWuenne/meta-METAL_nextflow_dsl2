# ## Test the nextflow workflow
# work_dir="$PWD"
# nextflow run main.nf -resume  \
# -with-report "../nextflow_reports/test_report.html" \
# -with-timeline "../nextflow_reports/test_timeline.html" \
# --meta_file "$work_dir/test_data/*.metal.txt" \
# --sumstat_files "./test_data/*mod.txt.gz" \
# --outdir "./test_out"


## Test with different input directories for sumstat files
work_dir="$PWD"
nextflow run main.nf -resume  \
-with-report "../nextflow_reports/test_report.html" \
-with-timeline "../nextflow_reports/test_timeline.html" \
--meta_file "$work_dir/test_data/*.metal.txt" \
--sumstat_files ["./test_data/meta_dir_1/*mod.txt.gz","./test_data/meta_dir_2/*mod.txt.gz"] \
--outdir "./test_out"