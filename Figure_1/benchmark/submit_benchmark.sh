#!/bin/bash

outdir=/nfs/research/marioni/akluzer/milo_analysis_2022/Figure_1/benchmark/data/output_benchmark

## Activate conda env
#conda activate milo_bm
source /hps/software/users/marioni/akluzer/miniconda3/etc/profile.d/conda.sh
conda activate milo_bm2

data_id=embryo

if [[ "$data_id" == "embryo" ]]
then 
pop_file=/nfs/research/marioni/akluzer/milo_analysis_2022/Figure_1/benchmark/data/pop_sample_2.txt
pops=$(cat $pop_file | while read pop; do echo $pop; done)
R_methods=$(for m in milo_red milo_tri; do echo $m; done)
batch_vec=$(for m in 0; do echo $m; done)
k=50
fi

output_dir=/nfs/research/marioni/akluzer/milo_analysis_2022/Figure_1/benchmark/data/output_benchmark

## Run
for pop in $pops
do
for pop_enr in $(seq 0.7 0.05 0.95)
do
for seed in $(seq 43 1 45)
do
for batch_sd in $batch_vec
do
for method in $R_methods
do
echo "Rscript /nfs/research/marioni/akluzer/milo_analysis_2022/Figure_1/benchmark/run_DA_R.r /nfs/research/marioni/akluzer/milo_analysis_2022/Figure_1/benchmark/data/${data_id}_data_bm.RDS $method $seed $pop --pop_enrichment $pop_enr --data_id $data_id --k $k --batchEffect_sd $batch_sd" | \
bsub -o ${outdir}/milo_bm_${data_id}_${seed}_${method}_${pop}_${batch_sd}.out -e ${outdir}/milo_bm_${data_id}_${seed}_${method}_${pop}_${batch_sd}.err -R"select[mem>50000] rusage[mem=50000]" -M50000 
#                             fi
done
done
done
done
done