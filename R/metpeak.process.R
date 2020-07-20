.metpeak.process <- function(IP,INPUT,batch_id,minimal_counts_in_fdr=10){
  
  #   using approximate Newton's method to call peaks per gene (the fastest method)
  Ng = unique(batch_id)
  nip <- ncol(IP)
  nin <- ncol(INPUT)
  INPUT_mean <- rowMeans(INPUT)
  IP_mean <- rowMeans(IP)
  if (nip > nin) {
    avg_input <- round(matrix(rep(INPUT_mean,nip-nin),ncol=nip-nin))
    INPUT <- cbind(INPUT,avg_input) 
  }
  else if (nip < nin){
    avg_ip <- matrix(rep(IP_mean,nin-nip),ncol=nin-nip)
    IP <- cbind(IP,avg_ip) 
  }
  # initialize the variables
  pvalues <- rep(1,nrow(IP))
  for (ii in Ng){
    
    print(ii)
    flag <- batch_id==ii
    ip=as.matrix(IP[flag,])
    input=as.matrix(INPUT[flag,])    
    # to avoid 0 read count in the input     
    bA <- max(apply(input,2,median)) + 1
    
    # to avoid ip reads is smaller than a bin of input
    # using approximate Newton written in R to do the peak calling
    res <-.betabinomial.hmm(ip,input+bA)
    cl <- .betabinomial.vit(ip,input+bA,res$trans,res$alpha, res$beta)
    clust_mean = res[[1]]/(res[[1]]+res[[2]])
    maxID = which.max(clust_mean)
    
    # Peak region should have mean ratio over 0.5, meaning IP reads should be more than Input reads in peak region
    if (clust_mean[maxID] >= 0.4){
      peak_loci <- which(cl$class==maxID)
      peak_loci_end <- length(peak_loci) # no peak condition included
      if (peak_loci_end > 0){
        peak_id <- c(0,which( (peak_loci[-1] - peak_loci[-peak_loci_end] ) > 1 ),peak_loci_end)
        for (jj in 1:(length(peak_id)-1)){
          jjpeak <- peak_loci[peak_id[jj]+1]:peak_loci[peak_id[jj+1]]
          res[[5]][jjpeak,maxID] <- median(res[[5]][jjpeak,maxID])  #res[[5]] = res$postprob
        }
      }  
      pvalues[flag] = 1 - res[[5]][,maxID]
    }
    else{
      pvalues[flag] = 1
    }
    
    
    
    #plot the result
    #     dot <- 1:nrow(ip); matplot(cbind(ip,input),col = 1:(ncol(ip)*2));points(dot[cl$class==1],rep(max(ip)/2,sum(cl$class==1))) # test part    
    #     anno = .get.gene.anno(ii,ANNOTATION,ANNOTATION_BATCH_ID)
    #     title(anno$gene)
    
    # cluster with biggest beta mean parameters assinged as peak cluster
    
    
  }
  log.fdr=log(p.adjust(pvalues,method='fdr'))
  
  # with significant number of reads only
  ID=which( (IP_mean+INPUT_mean) > minimal_counts_in_fdr) #should be vector not matrix
  log.fdr_sig=log(p.adjust(pvalues[ID], method = "fdr"))
  log.fdr[ID]=log.fdr_sig
  # fold enrichment
  log.fc=log(IP_mean/(INPUT_mean+1))
  # if total reads smaller than IP reads do not recognize a peak
  #   ID <- which( log.fc <= 0 )
  #   log.fdr[ID] <- 1
  #   pvalues[ID] <- 1
  
  # output result
  PW=list(log.p=log(pvalues),log.fdr=log.fdr,log.fc=log.fc)
  return(PW)
  
}
