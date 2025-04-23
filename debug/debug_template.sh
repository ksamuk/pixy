# template for debugging user submitted vcfs and pop files
# tbd: context specificity (e.g. ignore bedfile if not present)


pixy="python pixy/__main__.py"
data_folder="tests/main/data"

# vcf_path="${data_folder}/ag1000_pixy_invariant_test.vcf.gz"
# pop_path="${data_folder}/ag1000_populations_file.txt"

# $pixy --debug \
# --stats pi --vcf $vcf_path  \
# --populations $pop_path \
# --window_size 100000 \
# --n_cores 4 \
# --output_folder . \
# --output_prefix all_pi \
# --bypass_invariant_check

pixy="python pixy/__main__.py"
data_folder="debug"
bed_path="${data_folder}/BEDfile.txt"
vcf_path="${data_folder}/subsampled_10K.vcf.gz"
pop_path="${data_folder}/popfile.txt"

echo "$pixy --debug \
--stats dxy \
--vcf $vcf_path  \
--populations $pop_path \
--bed_file $bed_path \
--output_folder debug \
--output_prefix test" 


--debug --stats dxy --vcf debug/subsampled_10K.vcf.gz  --populations debug/popfile.txt --bed_file debug/BEDfile.txt --output_folder debug --output_prefix test

--debug --stats dxy --vcf debug/subsampled_10K.vcf.gz  --populations debug/popfile.txt --window_size 10000 --output_folder debug --output_prefix test
