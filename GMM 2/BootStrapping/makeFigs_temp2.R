library(beanplot)
source("parseData.R", local = TRUE)
source("calcPlotDefect_SD.R", local = TRUE)
dir.create(file.path("Output"), showWarnings = FALSE)
dir.create(file.path("PDF"), showWarnings = FALSE)
dir.create(file.path("PNG"), showWarnings = FALSE)

mchans = c("Ch1","Ch2","Ch3","Ch4","Area")
chans=c("LAMA1","MTCO1","VDAC1","NDUFB8","Area") # standard for CD7
# chans=c("LAMA1","VDAC1","MTCO1","NDUFB8","Area") # swapped MTCO1 and VDAC1 for confocal

fulldat =
  #"IMV.E02.P01.PMT.M3243AG.QIF.7TR.RAW.txt"
  #"IMV.E02.P01.PMT.M3243AG.QIF.TCF.RAW.txt"
  #"IMV.E01.P02.PMT.M3243AG.QIF.TGo.RAW.txt"
  "IMV.E02.P01.PMT.M3243AG.QIF.TGo.RAW.txt"
dat = read.delim(fulldat,stringsAsFactors=FALSE)
for(i in seq_along(mchans)){
  dat$Channel = gsub(mchans[i],chans[i],dat$Channel)
}

write.table(dat,"IMV.E02.P01.PMT.M3243AG.QIF.TGo.RAW_ren.txt",row.names=FALSE,quote=FALSE,sep="\t")

cord = c("NDUFB8","MTCO1","VDAC1")
chlabs = c("CI","CIV","OMM")
names(chlabs) = cord
mitochan = "VDAC1"
correctnpc = TRUE
updatechans = FALSE

dat = getData("IMV.E02.P01.PMT.M3243AG.QIF.TGo.RAW_ren.txt", cord,
              mitochan = mitochan,updatechans = updatechans, correctnpc = correctnpc)
dat$fn = gsub("_.0", "", dat$filename)
dat$pch = paste(dat$fn,dat$ch,sep="_")

# if not correcting, NPC need to be removed
#if(!correctnpc) dat = dat[grepl("OXPHOS",dat$filename),]

# grabbing and seperating crls and patients
sbj = sort(unique(dat$fn))
crl = grep("C._H", sbj, value = TRUE)
pts = grep("P.", sbj, value = TRUE)

# Function to determine data ranges to be used for graphs's axis
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

# scatter plots mito class with predictive interval
png("PNG/Fig3A.png",width=4500,height=6500,pointsize=60)
op = par(mfrow=c(3,2),mar = c(5.5,5.5,3,3))
i = 1
for(pat in pts[1:3]){
 ci = mitoplot(dat,pat,"NDUFB8",lab_inner=chlabs["NDUFB8"],xlim=gxylim,ylim=gxylim)
 labelpanel(paste0(tolower(as.roman(i)),")"))
 i = i+1
 civ = mitoplot(dat,pat,"MTCO1",lab_inner=chlabs["MTCO1"],xlim=gxylim,ylim=gxylim)
 labelpanel(paste0(tolower(as.roman(i)),")"))
 i = i+1
}
par(op)
dev.off()

png("PNG/Fig3B.png",width=4500,height=6500,pointsize=60)
op = par(mfrow=c(3,2),mar = c(5.5,5.5,3,3))
i = 1
for(pat in pts[4:6]){
  ci = mitoplot(dat,pat,"NDUFB8",lab_inner=chlabs["NDUFB8"],xlim=gxylim,ylim=gxylim)
  labelpanel(paste0(tolower(as.roman(i)),")"))
  i = i+1
  civ = mitoplot(dat,pat,"MTCO1",lab_inner=chlabs["MTCO1"],xlim=gxylim,ylim=gxylim)
  labelpanel(paste0(tolower(as.roman(i)),")"))
  i = i+1
}
par(op)
dev.off()

png("PNG/Fig4.png",width=15000,height=25000,pointsize=60)
op = par(mfrow=c(10,6),mar=c(5.5,5.5,3,3))
pats = pts
labs = LETTERS[0:length(pats)]
names(labs) = pats
for(pat in pts){
 ci = mitoplot(dat,pat,"NDUFB8",lab_inner=chlabs["NDUFB8"],xlim=gxylim,ylim=gxylim)
 labelpanel(labs[pat])
 civ = mitoplot(dat,pat,"MTCO1",lab_inner=chlabs["MTCO1"],xlim=gxylim,ylim=gxylim)
}
par(op)
dev.off()

getDeficiency = function(dat,pat,chan,bootstrap){
  res = mitoplot(dat,pat,chan,makeplot=FALSE,bootstrap=bootstrap)
  return(100*res$low/res$N)
}

Nboots = 2500
patboot = list()
targprots = cord[cord!=mitochan]

system.time(
for(pat in pats){
 bdf = data.frame(sample=1:Nboots)
 print(pat)
 for(ch in targprots) {print(ch); bdf[[ch]] = replicate(Nboots,getDeficiency(dat,pat,ch,TRUE))}
 bdf$sample = NULL
 patboot[[pat]] = bdf
 fpath = file.path("Output",paste0(pat,"_single_defect.csv"))
 write.table(patboot[[pat]],fpath,row.names=FALSE,quote=FALSE,sep=",",col.names = !file.exists(fpath), append = T)
}
)

base::saveRDS(patboot, "Output/patboot.Rds") # save bootstrap list
patboot1 = base::readRDS("Output/patboot.Rds") # load bootstrap list

 # Draw some plots, with dots
for(pat in pats){ 
 png(paste0("PNG/PlotsWithDots_",pat,".png"),width=9000,height=4500,pointsize=60)
 op = par(mfrow=c(1,2),mar=c(5.5,5.5,5.5,2.5))
  for(ch in targprots) mitoplot(dat,pat,ch,makeplot=TRUE,bootstrap=FALSE)
 par(op)
 dev.off()
}

png("PNG/StripBySection.png",width=3*3000/2,height=15000,pointsize=60)
op = par(mfrow=c(10,3),mar=c(6.5,5.5,5.5,2.5))
# Plot for each section
for(pat in pats){
 mlab = pat
 bdf = patboot[[pat]]
 beanplot(bdf,what=c(0,0,0,0),col=rgb(0,0,1,0.45),axes=TRUE,linelwd=3,border=NA,ylim=c(0,100),ylab="Deficient fibres (%)",xlab="",main=mlab,
   cex.axis=1.5,cex.lab=2.55,cex.main=2.55,las=2,wd=1.5,log="")
 abline(v=c(1.5),lty=2,lwd=4)
 text(x=c(1.0,2.0),y=100,c("CI", "CIV"),cex=1.5)
 stripchart(bdf,vertical=TRUE,pch=16,col=rgb(1,0,0,0.25),cex=0.35,method="jitter",jitter=0.3,add=TRUE)
}
par(op)
dev.off()

png("PNG/StripByProt.png",width=12000,height=8000,pointsize=60)
op = par(mfrow=c(2,1),mar=c(12.5,5.5,3,3))
# Plot for each block
for(ch in targprots){
  pdf = data.frame(sample=1:Nboots)
  for(pat in pats) pdf[[pat]] = patboot[[pat]][[ch]]
  pdf$sample=NULL
  mlab = paste0(ch," (",chlabs[ch],")")
  beanplot(pdf,what=c(0,0,0,0),col=rgb(0,0,1,0.45),axes=TRUE,linelwd=3,border=NA,ylim=c(0,100),ylab="Deficient fibres (%)",xlab="",main=mlab,
   cex.axis=1.5,cex.lab=2.55,cex.main=2.55,las=2,wd=1.5,log="")
  #text(x=c(1.75,5.05,7.65),y=100,c("class I", "class II", "class III"),cex=2.0)
  #abline(v=c(3.5,6.5),lty=2,lwd=4)
  stripchart(pdf,vertical=TRUE,pch=16,col=rgb(1,0,0,0.25),cex=1,method="jitter",jitter=0.3,add=TRUE) 
}
par(op)
dev.off()


mitop = dat[(dat$ch==mitochan)&(dat$type=="Mean intensity"),]

png("PNG/Histograms_VDAC1.png",width=6000,height=4500,pointsize=60)
op = par(mfrow=c(3,3),mar=c(5.5,5.5,5.5,2.5))
# Plot for each section
for(pat in pats){
 pclass = "m.3243A>G"
 #if(pat%in%c("P01","P02","P03")) pclass = "class I"
 #if(pat%in%c("P04","P05","P06")) pclass = "class II"
 #if(pat%in%c("P07","P08")) pclass = "class III"
 mlab = paste0(pat," (",pclass,")")
 mp = mitop$value[mitop$fn==pat]
 plot(density(log(mp)),main=mlab,xlim=range(log(mitop$value)),lwd=5,cex.axis=2.0,cex.lab=2.35,cex.main=2.55)
}
plot.new()
text(0.5,0.5,paste0("log(",mitochan,")"),cex=4.5)
par(op)
dev.off()

# environment backup if needed
save.image(file = "enviback.Rdata")
load("enviback.Rdata")

