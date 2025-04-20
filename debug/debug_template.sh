# template for debugging user submitted vcfs and pop files
# tbd: context specificity (e.g. ignore bedfile if not present)


pixy="python pixy/__main__.py"
data_folder="tests/main/data"

vcf_path="${data_folder}/ag1000_pixy_invariant_test.vcf.gz"
pop_path="${data_folder}/ag1000_populations_file.txt"

$pixy --debug \
--stats pi --vcf $vcf_path  \
--populations $pop_path \
--window_size 100000 \
--n_cores 4 \
--output_folder . \
--output_prefix all_pi \
--bypass_invariant_check
