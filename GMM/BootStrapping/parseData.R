library(data.table)

getData = function(fname,cord,mitochan="VDAC1",updatechans=TRUE,correctnpc=FALSE){
  dat = fread(fname,sep="\t",stringsAsFactors=FALSE,header=TRUE)
  colnames(dat) = tolower(colnames(dat))
  
  if(updatechans) {
    dat$channel = gsub("GRIM19","NDUFA13",dat$channel)
    dat$channel = gsub("Ch1","Laminin",dat$channel)
    dat$channel = gsub("Ch2","MTCO1",dat$channel)
    dat$channel = gsub("Ch3","VDAC1",dat$channel)
    dat$channel = gsub("Ch4","SDHA",dat$channel)
    dat$ch = substring(dat$channel,regexpr("\\_[^\\_]*$", dat$channel)+1,nchar(dat$channel))
  }else{
    dat$ch = as.character(dat$channel)
  }
  print("Available channels:")
  print(unique(dat$ch))
  dat = dat[dat$ch%in%cord,]
  dat$fn = gsub(" NPC","",dat$filename)
  dat$fn = gsub(" OXPHOS","",dat$fn)
  dat$pch = paste(dat$fn,dat$ch)
  
  if (correctnpc){
    npc = dat[grepl("NPC",dat$filename),]
    oxphos = dat[grepl("OXPHOS",dat$filename),]
    agg = aggregate(npc$value,by=list(npc$pch),mean,data=dat)
    lu = agg$x
    names(lu) = agg$Group.1
    
    oxphos$filename = gsub(" OXPHOS","",oxphos$filename)
    oxphos$value = pmax(1.0,oxphos$value - lu[oxphos$pch])
    oxphos$filename = oxphos$fn
    oxphos$fn = NULL
    oxphos$pch = NULL
    dat = oxphos
  }
  
  dat$type = "Mean intensity"
  dat$type[grepl("LOG_",dat$channel)] = "Log mean intensity"
  dat$type[grepl("MED_",dat$channel)] = "Median intensity"
  dat$type[grepl("R_",dat$channel)] = "Ratio mean intensity (VDAC1)"
  dat$type[grepl("R_MED_",dat$channel)] = "Ratio median intensity (VDAC1)"
  dat$type[grepl("R_LOG_",dat$channel)] = "Ratio log mean intensity (VDAC1)"
  dat$type[grepl("Z_",dat$channel)] = "z-score"
  dat$outlier_diff = "NODIFF"
  dat$regression_diff = "NODIFF"
  dat$z_diff = "NODIFF"
  dat$z = 0
  
  dat$chstr = dat$ch
  transform = log
  dat_r = dat[dat$type=="Mean intensity",]
  dat_r$type = paste("r (",mitochan,")",sep="")
  dat_theta = dat[dat$type=="Mean intensity",]
  dat_theta$type = paste("theta (",mitochan,")",sep="")
  
  for(pid in unique(dat$patrep_id)){
    for(ch in unique(dat$ch)){
      dt = dat[(dat$patrep_id==pid)&(dat$type=="Mean intensity"),]
      
      isch = as.character(dt$ch)==ch
      ismito = as.character(dt$ch)==mitochan
      prot = dt[isch,]
      mito = dt[ismito,]
      
      x = mito$value
      y = prot$value
      dat_r$value[(dat_r$patrep_id==pid)&(as.character(dat_r$ch)==ch)] = sqrt(x^2+y^2)
      dat_r$channel[(dat_r$patrep_id==pid)&(as.character(dat_r$ch)==ch)] = paste("RADIUS",ch,sep="_")
      dat_theta$value[(dat_theta$patrep_id==pid)&(as.character(dat_theta$ch)==ch)] = 360*atan(y/x)/(2*pi)
      dat_theta$channel[(dat_theta$patrep_id==pid)&(as.character(dat_theta$ch)==ch)] = paste("THETA",ch,sep="_")
    }
  }
  dat=rbind(dat,dat_r,dat_theta)
  dat = data.frame(dat, stringsAsFactors=FALSE)
  return(dat)
}