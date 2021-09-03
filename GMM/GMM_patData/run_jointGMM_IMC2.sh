#!/bin/sh

nohup Rscript bayes_jointGMM_IMC2.R OSCP &
nohup Rscript bayes_jointGMM_IMC2.R SDHA &
nohup Rscript bayes_jointGMM_IMC2.R NDUFB8 &
nohup Rscript bayes_jointGMM_IMC2.R GRIM19 &
nohup Rscript bayes_jointGMM_IMC2.R COX4+4L2 &
nohup Rscript bayes_jointGMM_IMC2.R UqCRC2 &
nohup Rscript bayes_jointGMM_IMC2.R MTCO1 &
