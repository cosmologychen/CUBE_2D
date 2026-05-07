lc_path='/mnt/18T/output_2D/lc_c10/'
run_dir='3000_1024_'
a_type=('runtime' 'runtime')
# b_type=('runtime' 'slide')
b_type=('onerun' 'one')




if [ ${a_type[1]} == 'slide' ]; then
    a_file="${lc_path}${run_dir}${a_type[0]}/lightcone_post_m1_xp.bin"
elif [ ${a_type[1]} == 'runtime' ]; then
    a_file="${lc_path}${run_dir}${a_type[0]}/runtime_lightcone_xp.bin"
elif [ ${a_type[1]} == 'one' ]; then
    a_file="${lc_path}${run_dir}${a_type[0]}/one_run_xp.bin"
    b_type=("one" "run")
fi
if [ ${b_type[1]} == 'slide' ]; then
    b_file="${lc_path}${run_dir}${b_type[0]}/lightcone_post_m1_xp.bin"
elif [ ${b_type[1]} == 'runtime' ]; then
    b_file="${lc_path}${run_dir}${b_type[0]}/runtime_lightcone_xp.bin" 
elif [ ${b_type[1]} == 'one' ]; then
    b_file="${lc_path}${run_dir}${b_type[0]}/one_run_xp.bin"
    b_type=("one" "run")
fi
output_file="${lc_path}comparsion/${a_type[0]}_${a_type[1]}_2_${b_type[0]}_${b_type[1]}.bin"
mkdir -p ${lc_path}comparsion/
# echo $a_file $b_file $output_file
# ls $a_file
# ls $b_file

make  && ./utilities/test_lc_geom.x $a_file $b_file $output_file