mchans = c("115In_CYB", "153Eu_SDHA", "160Gd_ND2", "161Dy_ATP8", "164Dy_ND4", 
"166Er_VDAC1", "168Er_COX4", "170Yb_ATPB", "172Yb_MTCO1", "174Yb_UqCRC2", 
"176Yb_Dystrophin")
chans=c("MTCYB","SDHA","MTND2","MTATP8","MTND4","VDAC1","COX4+4L2","MTATPB","MTCO1","UQCRC2","Dystrophin")

fulldat = "../mitocyto_merged_results.txt"
dat = read.delim(fulldat,stringsAsFactors=FALSE)
for(i in seq_along(mchans)){
  dat$Channel = gsub(mchans[i],chans[i],dat$Channel)
}

write.table(dat,"../mitocyto_merged_results_2.txt",row.names=FALSE,quote=FALSE,sep="\t")