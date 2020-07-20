

# comparing two experimental conditions
.report.peak<- function(TOTAL_PEAK_RESULT,PARAMETERS) {
  
  # get dir
  dir=paste(PARAMETERS$OUTPUT_DIR,PARAMETERS$EXPERIMENT_NAME,sep='/')
  
  write.table(TOTAL_PEAK_RESULT,file=paste(dir,"peak.xls",sep="/"), sep="\t",row.names =FALSE,quote = FALSE)
  temp=TOTAL_PEAK_RESULT[,1:12]
  names(temp)[1]=paste("#",names(temp)[1])
  write.table(temp,file=paste(dir,"peak.bed",sep="/"), sep="\t",row.names =FALSE,quote = FALSE)
  
}

# comparing two experimental conditions
report.peak<- function(TOTAL_PEAK_RESULT,PARAMETERS) {
  
  # get dir
  dir=paste(PARAMETERS$OUTPUT_DIR,PARAMETERS$EXPERIMENT_NAME,sep='/')
  
  write.table(TOTAL_PEAK_RESULT,file=paste(dir,"peak.xls",sep="/"), sep="\t",row.names =FALSE,quote = FALSE)
  temp=TOTAL_PEAK_RESULT[,1:12]
  names(temp)[1]=paste("#",names(temp)[1])
  write.table(temp,file=paste(dir,"peak.bed",sep="/"), sep="\t",row.names =FALSE,quote = FALSE)
  
}