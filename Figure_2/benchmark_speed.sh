#!/bin/bash

## Activate conda env
## source activate base
source /hps/software/users/marioni/akluzer/miniconda3/etc/profile.d/conda.sh
conda activate milo_bm

## Run embryo
echo "Rscript /nfs/research/marioni/akluzer/milo_analysis_2022/Figure_2/benchmark_speed.R" | bsub -q production -o /nfs/research/marioni/akluzer/milo_analysis_2022/Figure_2/speed_bm.out -e /nfs/research/marioni/akluzer/milo_analysis_2022/Figure_2/speed_bm.err -R"select[mem>70000] rusage[mem=70000]" -M70000
