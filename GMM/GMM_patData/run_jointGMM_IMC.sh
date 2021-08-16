#!/bin/sh

nohup Rscript bayes_jointGMM_IMC.R OSCP &
nohup Rscript bayes_jointGMM_IMC.R SDHA &
nohup Rscript bayes_jointGMM_IMC.R NDUFB8 &
nohup Rscript bayes_jointGMM_IMC.R GRIM19 &
nohup Rscript bayes_jointGMM_IMC.R COX4+4L2 &
nohup Rscript bayes_jointGMM_IMC.R UqCRC2 &
nohup Rscript bayes_jointGMM_IMC.R MTCO1 &
