library(rjags)
library(beanplot)
library(MASS)
source("../BootStrapping/parseData.R", local = TRUE)


for(fulldat in c("IMC.RAW.txt") ){
  # removes '.RAW.txt' from fulldat, stores as froot
  froot = gsub(".txt","",fulldat)
  
  # mchans = c("Ch1","Ch2","Ch3","Ch4","Area")
  # if(grepl(".TCF.", fulldat)){# attention confocal - swap channels
  #   chans=c("LAMA1","VDAC1","MTCO1","NDUFB8","Area") # swapped MTCO1 and VDAC1 for confocal
  # }else{
  #   chans=c("LAMA1","MTCO1","VDAC1","NDUFB8","Area") # standard for CD7
  # }
  # 
  # # names each element of chans with elements of mchans
  # names(chans) = mchans
  # 
  dat = read.delim(file.path("../BootStrapping",fulldat), stringsAsFactors=FALSE)
  # maps {Ch1,Ch2,Ch3,Ch4,Area} to {LAMA1,MTCO1,VDAC1,NDUFB8,Area}
  # dat$channel = chans[dat$Channel] 
  
  fulldat_ren = gsub(".RAW", ".RAW_ren", fulldat)
  
  write.table(dat, fulldat_ren, row.names=FALSE, 
              quote=FALSE, sep="\t")
  
  cord = c("NDUFA13", "NDUFB8", "SDHA", "UqCRC2", "COX4+4L2", "MTCO1",  
           "OSCP", "VDAC1") 
  # chlabs = c("CI","CIV","OMM") # chlabs - only used here. Not important?
  # names(chlabs) = cord # unsure
  mitochan = "VDAC1"
  correctnpc = TRUE
  updatechans = FALSE
  
  dat = getData( fulldat_ren, cord, mitochan = mitochan, updatechans = updatechans, 
                 correctnpc = correctnpc )
  
  dat$fn = gsub("_.0", "", dat$filename) # what is fn ??
  dat$pch = paste(dat$fn, dat$ch,sep="_") # what is pch ??
  
  # Get plot axis ranges
  lims = list()
  for(ch in cord){
    lims[[ch]] = quantile(log(dat$value[dat$channel==ch]), c(0.001,0.999),na.rm=TRUE)
  }
  
  # Merge different regions of the same section
  dat$fn = gsub("_R1","",dat$fn)
  dat$fn = gsub("_R2","",dat$fn)
  
  # grabbing and seperating ctrls and patients
  sbj = sort(unique(dat$fn))
  crl = grep("C._H", sbj, value = TRUE)
  pts = NULL
  pts = grep("P.", sbj, value = TRUE) # why define pts twice?
  
  #for(chan in c("NDUFB8","MTCO1")){
  for(chan in c("NDUFB8")){ # only run inference for NDUFB8
    
    # froot: data name 
    outroot_ctrl = paste(froot,"CTRL",chan,sep="__")
    # saves posterior draws in "Output" file
    posterior_ctrl_file = file.path("Output",paste0(outroot_ctrl,"__POSTERIOR.txt"))
    
    # dataframe for mean intensity of Control data 
    control = dat[(dat$fn%in%crl)&(dat$type=="Mean intensity"),]
    
    # identifies and log transforms the explanatory variable (VDAC1) and 
    # response variable (NDFUB8 or MTOC1)
    Xctrl = log(control$value[control$channel==mitochan])
    Yctrl = log(control$value[control$channel==chan])
    XY_ctrl = cbind( Xctrl, Yctrl )
    
    # if the output doesn't exist then infer control parameters
    if(!file.exists(posterior_ctrl_file)){
      
      # define prior parameters
      mu1_mean = c(0,0)
      mu2_mean = c(0,0)
      mu1_prec = 0.25*diag(2)
      mu2_prec = 0.25*diag(2)
      U_1 = solve(matrix(c(2,0,0,2), ncol=2, nrow=2, byrow=TRUE))
      n_1 = 2
      U_2 = solve(matrix(c(2,0,0,2), ncol=2, nrow=2, byrow=TRUE))
      n_2 = 2
      alpha_p = 2
      beta_p = 2
      pi = 0
      
      # Bayesian inference using JAGS
      # Assume a prior centred around the true values (assume good estimate from control data)
      
      Nctrl = nrow(XY_ctrl)
      # parameter list for RJAGS
      data_ctrl = list(Y=XY_ctrl, N=Nctrl, mu1_mean=mu1_mean, mu1_prec=mu1_prec, 
                       mu2_mean=mu2_mean, mu2_prec=mu2_prec, n_1=n_1, n_2=n_2,
                       U_1=U_1, U_2=U_2, alpha_p=alpha_p, beta_p=beta_p, pi=pi)
      
      data_ctrl_priorpred = data_ctrl # same parameters used for prior prediction RJAGS code
      data_ctrl_priorpred$Y = NULL # removes for prior prediction RJAGS 
      data_ctrl_priorpred$N = 0 # N: number of observed control points, removes for prior prediction
      
      # run the JAGS model for prior prediction and control parameter inference
      model_ctrl=jags.model(textConnection(modelstring), data=data_ctrl) # no initial vals given -> draw from prior
      
      model_ctrl_priorpred=jags.model(textConnection(modelstring), data=data_ctrl_priorpred) # no initial vals given -> draw from prior
      
      update(model_ctrl, n.iter=MCMCUpdates)
      
      output_ctrl=coda.samples(model=model_ctrl,variable.names=c("mu","tau","Y_syn","z","probdiff"),
                               n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
      
      output_ctrl_priorpred=coda.samples(model=model_ctrl_priorpred,variable.names=c("mu","tau","Y_syn", "z", "probdiff"),
                                         n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
      
      posterior_ctrl = as.data.frame(output_ctrl[[1]])
      prior_ctrl = as.data.frame(output_ctrl_priorpred[[1]])
      
      colnames(posterior_ctrl) = colnames(output_ctrl[[1]])
      colnames(prior_ctrl) = colnames(output_ctrl_priorpred[[1]])
      
      par(mfrow=c(2,3))
      plot(output_ctrl[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                          "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                          "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                          "probdiff", "Y_syn[1]", "Y_syn[2]")])
      par(mfrow=c(1,1))
      
      summ_ctrl = summary(output_ctrl)
      classifs_ctrl = summ_ctrl$statistics[grepl("z",rownames(summ_ctrl$statistics)),"Mean"]
      
      # print MCMC output
      #print(summ_ctrl)
      #plot(output_ctrl)
      #autocorr.plot(output_ctrl)
      #pairs(as.matrix(output_ctrl))
      #crosscorr.plot(output_ctrl)
      
      # prior and posterior prediction for control data
      predpsumm_ctrl=summary(output_ctrl_priorpred)
      ctrlroot = paste(froot,"CONTROL",chan,sep="__") 
      
      pdf(file.path("PDF",paste0(ctrlroot,".pdf")),width=14,height=7)
      priorpost(prior=prior_ctrl, posterior=posterior_ctrl, 
                data=data_ctrl, classifs=classifs_ctrl, ctrl=TRUE)
      title(paste(froot,"CTRL"), line = -1, outer = TRUE)
      dev.off()
      
      
      write.table(posterior_ctrl[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                                    "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                                    "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                                    "probdiff", "Y_syn[1]", "Y_syn[2]")],
                  posterior_ctrl_file,row.names=FALSE,quote=FALSE)
      
    }else{ # if file exists load previous data
      
      posterior_ctrl = read.delim(posterior_ctrl_file, sep=" ",stringsAsFactors=FALSE)
      
      colnames(posterior_ctrl) = c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                                   "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                                   "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                                   "probdiff", "Y_syn[1]", "Y_syn[2]")
      
      ctrlroot = paste(froot,"CONTROL",chan,sep="__") 
      
      pdf(file.path("PDF",paste0(ctrlroot,".pdf")), width=14, height=7)
      priorpost(ctrl_data=XY_ctrl, prior=prior_ctrl, posterior=posterior_ctrl, data=data_ctrl, ctrl=TRUE,
                class_posterior=rep(0,nrow(posterior_ctrl)), classifs=classifs_ctrl)
      title(paste(froot,"CTRL"), line = -1, outer = TRUE)
      dev.off()
      
    }
    
    ###
    ### prior specification for patient data
    ###
    
    n_1 = 6 # degrees of freedom
    # define the expected value of the patient prior (prec_pred) be the mean of the control
    # posterior
    prec_pred = matrix( colMeans(posterior_ctrl[,c('tau[1,1,1]', 'tau[1,2,1]', 'tau[1,2,1]','tau[2,2,1]')]),
                        nrow=2, ncol=2, byrow=TRUE)
    
    # increase the covariance between 'x' and 'y', keep variances the same
    # Sigma = solve(prec_pred)
    # delta = matrix(c(1,0.7,0.7,1), ncol=2, nrow=2, byrow=TRUE)
    # re-define the expectation of the prior
    # prec_pred = solve( Sigma + delta )
    # define prior parameter
    U_1 = prec_pred/n_1
    n_2 = 3
    U_2 = solve( matrix( c(2,0,0,2), nrow=2, ncol=2, byrow=TRUE) )/n_2
    
    mu1_mean = colMeans( posterior_ctrl[,c('mu[1,1]','mu[2,1]')])
    mu1_prec = solve( matrix( c(1,0,0,1), ncol=2, nrow=2, byrow=TRUE) )
    
    mu2_mean = c(0,0)
    mu2_prec = solve( matrix( c(5,0,0,5), ncol=2, nrow=2, byrow=TRUE) )
    
    alpha_p = 2
    beta_p = 2
    pi = 1
    
    for(pat in pts){ # loop through patients
      outroot = paste(froot,pat,chan,sep="__")
      patient = dat[(dat$fn==pat)&(dat$type=="Mean intensity"),] 
      posterior_file = file.path("Output",paste0(outroot,"__POSTERIOR.txt"))
      
      if(!file.exists(posterior_file)){ # regression for mitochondrial disease patients
        # Block off file from analysis
        file.create(posterior_file)
        
        op = par(mfrow=c(2,3), mar = c(5.5,5.5,3,3))
        
        Xpat = log(patient$value[patient$channel==mitochan])
        Ypat = log(patient$value[patient$channel==chan]) 
        XY_pat = cbind(Xpat, Ypat)
        
        # Bayesian inference using JAGS                                                                                                                                                  
        # Assume a prior centered around the estimate from control data
        Npat = nrow(XY_pat)
        
        data_pat = list(Y=XY_pat, N=Npat, mu1_mean=mu1_mean, mu1_prec=mu1_prec, 
                        mu2_mean=mu2_mean, mu2_prec=mu2_prec, n_1=n_1, n_2=n_2,
                        U_1=U_1, U_2=U_2, alpha_p=alpha_p, beta_p=beta_p, pi=pi)
        
        data_pat_priorpred = data_pat
        data_pat_priorpred$Y = NULL
        data_pat_priorpred$N = 0
        
        model_pat=jags.model(textConnection(modelstring), data=data_pat, n.chains=1) 
        
        model_pat_priorpred=jags.model(textConnection(modelstring), data=data_pat_priorpred) 
        update(model_pat,n.iter=MCMCUpdates)
        
        converge_pat = coda.samples(model=model_pat,variable.names=c("mu","tau","Y_syn","z","probdiff"),
                                    n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
        
        output_pat = coda.samples(model=model_pat,variable.names=c("mu", "tau","Y_syn","z","probdiff"),
                                  n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
        
        output_pat_priorpred = coda.samples(model=model_pat_priorpred,
                                            variable.names=c("mu", "tau","Y_syn","z","probdiff"),
                                            n.iter=MCMCUpdates_Report,thin=MCMCUpdates_Thin)
        
        plot(output_pat[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                           "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                           "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                           "probdiff", "Y_syn[1]", "Y_syn[2]")] )
        
        posterior_pat = as.data.frame(output_pat[[1]])
        prior_pat = as.data.frame(output_pat_priorpred[[1]])
        
        class_posterior_pat = posterior_pat[, grepl('z', colnames(posterior_pat))]
        colnames(posterior_pat) = colnames(output_pat[[1]])
        colnames(prior_pat) = colnames(output_pat_priorpred[[1]])
        
        summ_pat = summary(output_pat)
        #classifs_pat = summ_pat$statistics[grepl("z",rownames(summ_pat$statistics)),"Mean"]
        classifs_pat = colMeans(class_posterior_pat)
        
        #print(summ_pat)
        #plot(converge_pat)
        #autocorr.plot(converge_pat)
        #pairs(as.matrix(converge_pat))
        #crosscorr.plot(converge_pat)
        
        #predpsumm_pat=summary(output_pat_priorpred)
        pdf(file.path("PDF",paste0(outroot,".pdf")),width=14,height=8.5)
        priorpost(prior=prior_pat, posterior=posterior_pat, 
                  data=data_pat, classifs=classifs_pat)
        title(paste(froot,pat), line = -1, outer = TRUE)
        dev.off()
        write.table(as.numeric(classifs_pat),file.path("Output",paste0(outroot,"__CLASS.txt")),row.names=FALSE,quote=FALSE,col.names=FALSE)
        write.table(posterior_pat[,c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
                                     "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
                                     "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
                                     "probdiff", "Y_syn[1]", "Y_syn[2]")],posterior_file,row.names=FALSE,quote=FALSE)
      }else{ # if file exists load previous data
        
        # class_pat_file = file.path("Output", paste0(outroot, "__CLASS.txt"))
        # 
        # posterior_pat = read.delim(posterior_file,sep=" ",stringsAsFactors=FALSE)
        # class_pat = read.delim(class_pat_file, sep="\n", header=FALSE)
        # 
        # colnames(posterior_pat) = c("mu[1,1]","mu[1,2]","mu[2,1]","mu[2,2]",
        #                             "tau[1,1,1]","tau[1,2,1]","tau[2,1,1]","tau[2,2,1]",
        #                             "tau[1,1,2]","tau[1,2,2]","tau[2,1,2]","tau[2,2,2]",
        #                             "probdiff", "Y_syn[1]", "Y_syn[2]")
        # 
        # outroot = paste(froot, pat, chan,sep="__") 
        # 
        # pdf(file.path("PDF",paste0(outroot,".pdf")), width=14, height=7)
        # priorpost(prior=prior_pat, posterior=posterior_pat,data=data_pat, 
        #           class_posterior=class_pat[,1], classifs=classifs_pat)
        # title(paste(froot, pat), line = -1, outer = TRUE)
        # dev.off()
      } 
      
    } # pats
  } # chans
} # fulldats
