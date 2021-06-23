library(beanplot)
source("parseData.R", local = TRUE)
source("calcPlotDefect.R", local = TRUE)
dir.create(file.path("Output"), showWarnings = FALSE)
dir.create(file.path("PDF"), showWarnings = FALSE)
dir.create(file.path("PNG"), showWarnings = FALSE)

fulldats = c(
  "IMV.E02.P01.PMT.M3243AG.QIF.7TR.RAW.txt",
  "IMV.E01.P02.PMT.M3243AG.QIF.TGo.RAW.txt",
  "IMV.E02.P01.PMT.M3243AG.QIF.TGo.RAW.txt",
  "IMV.E02.P01.PMT.M3243AG.QIF.TCF.RAW.txt"
)

for(fulldat in fulldats){
print(fulldat)

froot = gsub(".RAW.txt","",fulldat)

Nboots = 250

mchans = c("Ch1","Ch2","Ch3","Ch4","Area")
if(grepl(".TCF.",fulldat)){# attention confocal - swap channels
  chans=c("LAMA1","VDAC1","MTCO1","NDUFB8","Area") # swapped MTCO1 and VDAC1 for confocal
}else{
  chans=c("LAMA1","MTCO1","VDAC1","NDUFB8","Area") # standard for CD7
}

names(chans) = mchans
pint = c(1:8)

dat = read.delim(fulldat,stringsAsFactors=FALSE)
dat$Channel = chans[dat$Channel]

write.table(dat,gsub(".RAW", ".RAW_ren", fulldat),row.names=FALSE,quote=FALSE,sep="\t")

cord = c("NDUFB8","MTCO1","VDAC1")
chlabs = c("CI","CIV","OMM")
names(chlabs) = cord
mitochan = "VDAC1"
correctnpc = TRUE
updatechans = FALSE

dat = getData(gsub(".RAW", ".RAW_ren", fulldat), cord,
              mitochan = mitochan,updatechans = updatechans, correctnpc = correctnpc)
dat$fn = gsub("_.0", "", dat$filename)
dat$pch = paste(dat$fn,dat$ch,sep="_")

# if not correcting, NPC need to be removed
#if(!correctnpc) dat = dat[grepl("OXPHOS",dat$filename),]

# grabbing and seperating ctrls and patients
sbj = sort(unique(dat$fn))
crl = grep("C._H", sbj, value = TRUE)
pts = grep("P.", sbj, value = TRUE)

# Function to determine data ranges to be used for graph axes
getRange = function(x, transform = log, quants = c(0.0,1.0)){ # function will transform x according to arguments (i.e. log)
  vals = transform(x) # it will log transform x to a vals vector
  vals = vals[!is.na(vals)&is.finite(vals)] # it will remove NaN and Inf values (not good for range calculations)
  return(as.numeric(quantile(vals,probs=quants))) #it will return just the min and max values via range() function
} #very important, close the function loop and activates return

xlim = getRange(dat$value[dat$ch==mitochan])
xlim
ylim_ci = getRange(dat$value[dat$ch=="NDUFB8"])
ylim_ci
ylim_civ = getRange(dat$value[dat$ch=="MTCO1"])
ylim_civ
gxylim = getRange(dat$value)
gxylim

# Demonstrating varying predictive interval definition:
pisdres = list()
for(pat in pts){
 print(pat)
 #pisds = seq(0,20,0.2)
 #defci = rep(0,length(pisds))
 #defciv = rep(0,length(pisds))
 #for(i in seq_along(pisds)){
 # pisd = pisds[i]
 # ci = mitoplot(dat,pat,"NDUFB8",lab_inner=chlabs["NDUFB8"],xlim=gxylim,ylim=gxylim,makeplot=FALSE,pisd=pisd)
 # civ = mitoplot(dat,pat,"MTCO1",lab_inner=chlabs["MTCO1"],xlim=gxylim,ylim=gxylim,makeplot=FALSE,pisd=pisd)
 # defci[i] = ci$low/ci$N
 # defciv[i] = civ$low/civ$N
 #}

 # Rate of change
 #diffci = vector(length=length(pisds))
 #diffciv = vector(length=length(pisds))
 #for(i in 2:(length(pisds)-1)) {
 # diffci[i] = 100*(defci[i+1]-defci[i-1])/(pisds[i+1]-pisds[i-1])
 # diffciv[i] = 100*(defciv[i+1]-defciv[i-1])/(pisds[i+1]-pisds[i-1])
 ##}
 # Find last time that the rate of change of proportion dips under 5% per SD
 #cutci = Position(isFALSE,diffci>-5,right=TRUE)
 #sdci = (pisds[cutci]+pisds[cutci+1])/2
 #cutciv = Position(isFALSE,diffciv>-5,right=TRUE)
 #sdciv = (pisds[cutciv]+pisds[cutciv+1])/2
 #pisdres[[pat]] = list(pisds=pisds,defci=defci, defciv=defciv,diffci=diffci,diffciv=diffciv,cutci=cutci,sdci=sdci,cutciv=cutciv,sdciv=sdciv)
 pisdres[[pat]] = list(sdci=qnorm(0.975),sdciv=qnorm(0.975))
}

pdf(paste0("PDF/VarPred_",froot,".pdf"),width=14,height=7)
for(pat in pts){
 pisds = pisdres[[pat]]$pisds
 defci = pisdres[[pat]]$defci
 defciv = pisdres[[pat]]$defciv
 diffci = pisdres[[pat]]$diffci
 diffciv = pisdres[[pat]]$diffciv
 cutci = pisdres[[pat]]$cutci
 cutciv = pisdres[[pat]]$cutciv
 sdci = pisdres[[pat]]$sdci
 sdciv = pisdres[[pat]]$sdciv

 op = par(mfrow=c(1,2),mar = c(5.5,5.5,3,3))
 plot(pisds,100*defci,type="l",lwd=3,col="blue",xlab = "Predictive interval range (SDs)",ylab = "Proportion of deficient fibres (%)",cex.axis=1.55,cex.lab=1.45,main=pat,ylim=c(0,100),xlim=c(0,10))
 points(pisds,100*defciv,type="l",lwd=3,col="red")
 abline(v=qnorm(0.975),lty=2,lwd=2,col="black")
 abline(v=c(sdci,sdciv),col=c("blue","red"),lwd=2,lty=2)
 legend("topright",c("CI","CIV"),col=c("blue","red"),lwd=3,cex=2)

 plot(pisds,diffci,type="l",lwd=3,col="blue",xlab = "Predictive interval range (SDs)",ylab = "Rate of change of proportion of deficient fibres (%/SD)",cex.axis=1.55,cex.lab=1.45,main=pat,ylim=c(-50,0),xlim=c(0,10))
 points(pisds,diffciv,type="l",lwd=3,col="red")
 abline(v=qnorm(0.975),lty=2,lwd=2,col="black")
 abline(v=c(sdci,sdciv),col=c("blue","red"),lwd=2,lty=2)
 legend("bottomright",c("CI","CIV"),col=c("blue","red"),lwd=3,cex=2)
 par(op)
}
dev.off()


getDeficiency = function(dat,pat,chan,bootstrap,pisd=qnorm(0.975)){
  res = mitoplot(dat,pat,chan,makeplot=FALSE,bootstrap=bootstrap,pisd=pisd)
  return(100*res$low/res$N)
}

patboot = list()
targprots = cord[cord!=mitochan]

dobootstrap = FALSE
if(dobootstrap){
 system.time(
 for(pat in pts){
  bdf = data.frame(sample=1:Nboots)
  print(pat)
  for(ch in targprots) {print(ch); bdf[[ch]] = replicate(Nboots,getDeficiency(dat,pat,ch,TRUE,pisd=ifelse(ch=="NDUFB8",pisdres[[pat]]$sdci,pisdres[[pat]]$sdciv)))}
  bdf$sample = NULL
  patboot[[pat]] = bdf
  fpath = file.path("Output",paste0(pat,"_single_defect.csv"))
  #write.table(patboot[[pat]],fpath,row.names=FALSE,quote=FALSE,sep=",",col.names = !file.exists(fpath), append = T)
  write.table(patboot[[pat]],fpath,row.names=FALSE,quote=FALSE,sep=",")
 }
 )
}

# Load in from .csv files
patboot = list()
for(pat in pts){
  fpath = file.path("Output",paste0(pat,"_single_defect.csv"))
  patboot[[pat]] = read.csv(fpath,stringsAsFactors=FALSE)
}

 # Draw some plots, with dots
for(pat in pts){
 png(paste0("PNG/PlotsWithDots_",pat,".png"),width=9000/4,height=4500/4,pointsize=60/3)
 op = par(mfrow=c(1,2),mar=c(5.5,5.5,5.5,2.5))
  for(ch in targprots) mitoplot(dat,pat,ch,makeplot=TRUE,bootstrap=FALSE,pisd=ifelse(ch=="NDUFB8",pisdres[[pat]]$sdci,pisdres[[pat]]$sdciv))
 par(op)
 dev.off()
}

png(paste0("PNG/StripBySection_",froot,".png"),width=1500,height=1000,pointsize=12)
op = par(mfrow=c(2,3),mar=c(6.5,5.5,5.5,2.5))
# Plot for each section
for(pat in pts){
 mlab = pat
 bdf = patboot[[pat]]
 beanplot(bdf,what=c(0,0,0,0),col=rgb(0,0,1,0.45),axes=TRUE,linelwd=3,border=NA,ylim=c(0,100),ylab="Deficient fibres (%)",xlab="",main=mlab,
   cex.axis=1.5,cex.lab=2.55,cex.main=2.55,las=2,wd=1.5,log="")
 abline(v=c(1.5),lty=2,lwd=4)
 text(x=c(1.0,2.0),y=100,c("CI", "CIV"),cex=1.5)
 stripchart(bdf,vertical=TRUE,pch=16,col=rgb(1,0,0,0.25),cex=1,method="jitter",jitter=0.3,add=TRUE)
}
par(op)
dev.off()

#png(paste0("PNG/StripByProt_",froot,".png"),width=1500,height=1000,pointsize=12)
png(paste0("PNG/StripByProt_",froot,"%02d.png"),width=2000,height=1500,pointsize=24)
#op = par(mfrow = c(2,1),mar=c(12.5,5.5,3,3))
op = par(mar=c(10.5,5.5,5.5,2.5))
# Plot for each block
for(ch in targprots){
  pdf = data.frame(sample=1:Nboots)
  for(pat in pts) pdf[[pat]] = patboot[[pat]][[ch]]
  pdf$sample=NULL
  todrop = colnames(pdf)[grepl("_R2",colnames(pdf))]
  for(col in todrop) pdf[[col]]=NULL
  colnames(pdf) = gsub("P1_","",colnames(pdf))
  colnames(pdf) = gsub("P2_","",colnames(pdf))
  colnames(pdf) = gsub("_R1","",colnames(pdf))
  blockids = substr(colnames(pdf),1,nchar(colnames(pdf))-3)
  tissues = substr(colnames(pdf),1,2)
  blocks = unique(blockids)

  pdf = pdf/100.0
  
  mlab = paste0(froot," ",ch," (",chlabs[ch],")")
  beanplot(pdf,what=c(0,0,0,0),col=rgb(0,0,1,0.45),axes=FALSE,linelwd=3,border=NA,ylim=c(0.0,1.0),ylab="Proportion of deficient fibres",xlab="",main=mlab,
   cex.axis=1.5,cex.lab=2.55,cex.main=2.0,las=2,wd=1.5,log="")
  #text(x=c(1.75,5.05,7.65),y=100,c("class I", "class II", "class III"),cex=2.0)
  #abline(v=c(3.5,6.5),lty=2,lwd=4)
  stripchart(pdf,vertical=TRUE,pch=16,col=rgb(0,0,0,0.25),cex=0.75,method="jitter",jitter=0.3,add=TRUE)
  axis(1,at=seq(2,length(blocks)*3,by=3),labels=blocks,las=2,cex.axis=2.0)
  axis(2,cex.axis=2.0)
  abline(v=seq(3.5,length(blocks)*3,by=3),lty=2)
  
}
par(op)
dev.off()

mitop = dat[(dat$ch==mitochan)&(dat$type=="Mean intensity"),]

png(paste0("PNG/Histograms_VDAC1_",froot,".png"),width=6000,height=4500,pointsize=60)
op = par(mfrow=c(3,3),mar=c(5.5,5.5,5.5,2.5))
# Plot for each section
for(pat in pts){
 pclass = "m.3243A>G"
 mlab = paste0(pat," (",pclass,")")
 mp = mitop$value[mitop$fn==pat]
 plot(density(log(mp)),main=mlab,xlim=range(log(mitop$value)),lwd=5,cex.axis=2.0,cex.lab=2.35,cex.main=2.55)
}
plot.new()
text(0.5,0.5,paste0("log(",mitochan,")"),cex=4.5)
par(op)
dev.off()

}
