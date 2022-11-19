#!/bin/bash
work_path=$(dirname $0)
name=$(basename $work_path)
# echo `date +%Y%m%d%H%M%S`

p=$1

infer_file=$2
exp=$name

python -u tools/main.py \
    --load_best \
    --infer \
    --infer_file $infer_file \
    --p_name $p\
    --out_dir $work_path \
    --exp_name $exp\
    ${@:6}| tee $work_path/out/log.txt
