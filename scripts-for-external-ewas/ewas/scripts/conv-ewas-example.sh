#!/bin/bash


## The command in the Snakefile
# rule conv_ewas:
#     input:
#         script = LOCAL_DIR + "/" + "scripts/conventional-ewas.R", 
#         phen = WORK_DIR + "/" + "data/cleaned-pheno-data.tsv",
#         meth = WORK_DIR + "/" + "data/clean-meth.RData",
#         svs = WORK_DIR + "/" + "data/svs/ad-svs.tsv",
#         covars = WORK_DIR + "/" + "data/covars-{model}.txt",
#     output: 
#         WORK_DIR + "/" + "results/ewas/ewaff-res-{model}.tsv", 
#     shell:
#         """echo $HOSTNAME; Rscript {input.script} \
#                                         '{input.phen}' \
#                                         '{input.meth}' \
#                                         '{input.svs}' \
#                                         '{input.covars}' \
#                                         '{output}' """   


LOCAL_DIR="."
WORK_DIR="."

## arguments
script="${LOCAL_DIR}/scripts/conventional-ewas.R"
phen="${WORK_DIR}/data/cleaned-pheno-data.tsv"
meth="${WORK_DIR}/data/clean-meth.RData"
svs="${WORK_DIR}/data/svs/ad-svs.tsv"

## cell count adjusted EWAS specific arguments
covars="${WORK_DIR}/data/covars-cc.txt"
output="${WORK_DIR}/results/ewas/ewaff-res-cc.tsv"

## cell count adjusted EWAS script
Rscript ${script} "${phen}" "${meth}" "${svs}" "${covars}" "${output}"


## cell count UNadjusted EWAS specific arguments
covars="${WORK_DIR}/data/covars-no-cc.txt"
output="${WORK_DIR}/results/ewas/ewaff-res-no-cc.tsv"

## cell count UNadjusted EWAS script
Rscript ${script} "${phen}" "${meth}" "${svs}" "${covars}" "${output}"

