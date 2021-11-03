labelled_data = read.csv("CNR_labelled_IMC.csv", stringsAsFactors=F)

imc_chan = imc_chan = c('SDHA','OSCP', 'GRIM19', 'MTCO1', 'NDUFB8', 'COX4+4L2', 'UqCRC2')
pts = c("P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10")

tt = read.csv("./Output/IMC_joint2/GRIM19__P01__CLASS.txt", sep=" ")
dim(tt)

labelled_data[labelled_data[,"GRIM19_up"] || labelled_data[,"GRIM19_down"], "cell_id"]

class_comp = function(labelled_data){
  output = vector("list", length=length(imc_chan))
  for(chan in imc_chan){
    chan_mat = matrix(NA, nrow=length(pts), ncol=3,
                                   dimnames=list(pts, c("like ctrl", "not like ctrl", "unsure")) )
    for(pat in pts){
      bayes_class = read.csv(paste0("./Output/IMC_joint2/", chan,"__",pat,"__CLASS.txt"), sep=" ")
      colnames(bayes_class) = c("cell_id", "classif")
      cell_ids = bayes_class[,"cell_id"]
      classifs = bayes_class[order(bayes_class["cell_id"]),]
      nfib = nrow(bayes_class)
      
      labelled_subset = labelled_data[labelled_data$cell_id %in% cell_ids, ]                       
      like_ctrl = sum(bayes_class[classifs<0.05, "cell_id"] %in% labelled_subset[labelled_subset[,paste0(chan,"_norm")], "cell_id"] )
      notlike_ctrl = sum( sum(bayes_class[classifs>0.95, "cell_id"] %in% labelled_subset[labelled_subset[,paste0(chan,"_up")] || labelled_subset[,paste(chan,"_down")], "cell_id"] ) )
      dunno_ctrl = nfib - like_ctrl - notlike_ctrl
      
      chan_mat[i,] = c(like_ctrl, notlike_ctrl, dunno_ctrl)
    }
    output[[paste(chan)]] = chan_mat
  }
  output
}

class_comp(labelled_data)



