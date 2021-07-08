## Create new tests files containing chromosome number and name for plotting for each test data set

## DGI
zcat DGI_three_regions.txt.gz | gawk -F " " '{print $1":"$2":A:G",$0}' > DGI_three_regions.mod.txt
sed -i 's/CHR:POS:A:G/new_ID/g' DGI_three_regions.mod.txt
gzip DGI_three_regions.mod.txt

## MAGIC FUSION
zcat MAGIC_FUSION_Results.txt.gz | gawk -F " " '{print $1":"$2":A:G",$0}' > MAGIC_FUSION_Results.mod.txt
sed -i 's/CHR:POS:A:G/new_ID/g' MAGIC_FUSION_Results.mod.txt
gzip MAGIC_FUSION_Results.mod.txt

## magic_SARDINIA
zcat magic_SARDINIA.tbl.gz | gawk -F " " '{print $2":"$3":A:G",$0}' > magic_SARDINIA.mod.tbl
sed -i 's/CHR:POS:A:G/new_ID/g' magic_SARDINIA.mod.tbl
gzip magic_SARDINIA.mod.tbl