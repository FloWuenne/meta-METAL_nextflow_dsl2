# ## Test the nextflow workflow
work_dir="$PWD"
nextflow run main.nf -resume  \
-with-report "../nextflow_reports/test_report.html" \
-with-timeline "../nextflow_reports/test_timeline.html" \
--meta_file "./test_data/*metal.txt"





--phenoFile "$work_dir/test_data/input/pheno*.txt" \
--phenoCol "y_binary" \
--covarColList "x1,x2" \
--bgen_prefix "genotype_100markers.chr" \
--bgen_suffix ".bgen" \
--bgen_path "$work_dir/test_data/input" \
--sampleFile "$work_dir/test_data/input/samplefile_test_input.txt" \
--outdir "../saige_test_out" \
--gwas_cat "../gwascat.csv"