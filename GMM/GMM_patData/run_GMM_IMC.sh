#!/bin/bash

nohup Rscript bayes_GMM_IMC.R NDUFB8 & 
nohup Rscript bayes_GMM_IMC.R GRIM19 & 
nohup Rscript bayes_GMM_IMC.R SDHA & 
nohup Rscript bayes_GMM_IMC.R OSCP & 
nohup Rscript bayes_GMM_IMC.R MTCO1 &  
nohup Rscript bayes_GMM_IMC.R COX4+4L2 & 
nohup Rscript bayes_GMM_IMC.R UqCRc2 & 
