#!/bin/bash

pop_file=/nfs/research/marioni/akluzer/milo_analysis_2022/Figure_1/benchmark/data/pop_sample_2.txt

outdir=/nfs/research/marioni/akluzer/milo_analysis_2022/Figure_1/benchmark/data

## Activate conda env
## source activate base
source /hps/software/users/marioni/akluzer/miniconda3/etc/profile.d/conda.sh
conda activate milo_bm

## Run embryo
cat $pop_file | \
while read pop
    do
    for pop_enr in $(seq 0.7 0.05 0.95)
        do
        for seed in $(seq 43 1 45)
            do
            echo "Rscript ./make_bm_data.R --pop_enrichment $pop_enr /nfs/research/marioni/akluzer/milo_analysis_2022/Figure_1/benchmark/embryo_data_bm.RDS $seed $pop" | bsub -o ${outdir}/milo_make_bm_data_${seed}.out -e ${outdir}/milo_make_bm_data_${seed}.err -R"select[mem>8000] rusage[mem=8000]" -M8000
            done
        done
    done