### Run DA methods in R ###

suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(SingleCellExperiment)
  library(scran)
})

source('/nfs/research/marioni/akluzer/milo_analysis_2022/Figure_1/benchmark/benchmark_utils.R')

parser <- ArgumentParser()
parser$add_argument("data_RDS", type="character",
                    help = "path to RDS storing SingleCellExperiment object")
parser$add_argument("method", type="character",
                    help = "DA method to use")
parser$add_argument("batch_seed", type="integer",
                    help = "Seed 4 batch effect")
parser$add_argument("population", type="character",
                    help = "Which cell type is DA?")
parser$add_argument("--pop_enrichment", type="double", default=0.7,
                    help = "Max condition probability in DA population")
parser$add_argument("--batchEffect_sd", type="double", default=0,
                    help = "SD of batch effect")
parser$add_argument("--k", type="integer", default=50,
                    help = "K parameter")
parser$add_argument("--data_id", type="character", default="embryo",
                    help = "ID for the dataset used")
parser$add_argument("--MNN_correct", type="character", default="no",
                    help = "should I use MNN corrected version")
parser$add_argument("--capped_enr", type="character", default="no",
                    help = "should I use the capped enrichment version")
parser$add_argument("--logit_param", type="character", default="no",
                    help = "if capped_enr != 'no', the a logit parameter")
args <- parser$parse_args()

seed <- args$batch_seed
data_path <- args$data_RDS
k <- args$k
pop <- args$population
pop_enr <- args$pop_enrichment
be_sd <- args$batchEffect_sd
DA_method <- args$method
data_id <- args$data_id
mnn_correct <- args$MNN_correct
capped <- args$capped
a_logit <- args$logit_param

## Load data
print("Loading dataset...")
sce <- readRDS(data_path)

## Load coldata and PCA
if (str_detect(pop, "/")) {
  pop <- str_replace_all(pop, "/", "_")
}
outdir <- '/nfs/research/marioni/akluzer/milo_analysis_2022/Figure_1/benchmark/data/'
outprefix <- str_c("benchmark_", data_id, "_pop_", pop, '_enr', pop_enr, "_seed", seed)

if (capped!="no") {
  outprefix <- str_c("benchmark_", data_id, "_pop_", pop, '_enr', pop_enr, "_seed", seed, "_capped_logit", a_logit)
}
coldata <- read_csv(paste0(outdir, outprefix, ".coldata.csv")) %>% column_to_rownames()
if (mnn_correct=="no") {
  X_pca <-read_csv(str_c(outdir, str_remove(outprefix, "_capped_logit.."), "_batchEffect", be_sd, ".pca.csv")) %>% column_to_rownames()
} else {
  X_pca <-read_csv(str_c(outdir, str_remove(outprefix, "_capped_logit.."), "_batchEffect", be_sd, ".MNNcorrected.pca.csv")) %>% column_to_rownames()
}

## Find DA probability x cell
tol_dataset <- list(embryo=1)

bm_params = list(
  milo_red = list(k=k),
  milo_tri = list(k=k))

## Run DA method
out <- runDA(sce, X_pca, coldata = coldata, method = DA_method, params = bm_params)

## Save output 
print("saving output")
bm_outdir <-'/nfs/research/marioni/akluzer/milo_analysis_2022/Figure_1/benchmark/data/output_benchmark/'
if (mnn_correct=="no") {
  write_csv(out, str_c(bm_outdir, outprefix, "_batchEffect", be_sd, ".DAresults.", DA_method, ".csv"))
} else {
  write_csv(out, str_c(bm_outdir, outprefix, "_batchEffect", be_sd, ".MNNcorrected.DAresults.", DA_method, ".csv"))
} 
