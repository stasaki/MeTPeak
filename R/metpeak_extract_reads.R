metpeak_extract_reads <- function(
  IP_BAM,INPUT_BAM,
  GENOME = NA,
  UCSC_TABLE_NAME = "knownGene",
  GENE_ANNO_GTF = NA,
  TXDB = NA,
  OUTPUT_DIR=NA,
  EXPERIMENT_NAME="MeTPeak_output",
  WINDOW_WIDTH=50,
  SLIDING_STEP=50,
  FRAGMENT_LENGTH=100,
  READ_LENGTH=NA,
  MINIMAL_PEAK_LENGTH=FRAGMENT_LENGTH/2,
  PEAK_CUTOFF_PVALUE=NA,
  PEAK_CUTOFF_FDR=0.05,
  FOLD_ENRICHMENT=1,
  MINIMAL_MAPQ=30,
  REMOVE_LOCAL_TAG_ANOMALITIES=TRUE,
  POISSON_MEAN_RATIO=1,
  TESTING_MODE=NA,
  SAVE_INTERMEDIATE=TRUE,
  METHOD=0
  
){

  # Wrap parameters ##################################################
  PARAMETERS=list();
  PARAMETERS$GENE_ANNO_GTF=GENE_ANNO_GTF
  PARAMETERS$IP_BAM=IP_BAM
  PARAMETERS$INPUT_BAM=INPUT_BAM
  PARAMETERS$GENOME = GENOME
  PARAMETERS$UCSC_TABLE_NAME=UCSC_TABLE_NAME
  UCSC_TABLE_NAME = UCSC_TABLE_NAME
  PARAMETERS$FRAGMENT_LENGTH=FRAGMENT_LENGTH
  PARAMETERS$READ_LENGTH=READ_LENGTH
  PARAMETERS$WINDOW_WIDTH=WINDOW_WIDTH
  PARAMETERS$SLIDING_STEP=SLIDING_STEP
  PARAMETERS$MINIMAL_PEAK_LENGTH=MINIMAL_PEAK_LENGTH
  PARAMETERS$PEAK_CUTOFF_PVALUE=PEAK_CUTOFF_PVALUE
  PARAMETERS$PEAK_CUTOFF_FDR=PEAK_CUTOFF_FDR
  PARAMETERS$FOLD_ENRICHMENT=FOLD_ENRICHMENT
  PARAMETERS$MINIMAL_MAPQ=MINIMAL_MAPQ
  PARAMETERS$REMOVE_LOCAL_TAG_ANOMALITIES=REMOVE_LOCAL_TAG_ANOMALITIES
  PARAMETERS$POISSON_MEAN_RATIO=POISSON_MEAN_RATIO
  PARAMETERS$OUTPUT_DIR=OUTPUT_DIR
  PARAMETERS$EXPERIMENT_NAME=EXPERIMENT_NAME
  PARAMETERS$TESTING_MODE=TESTING_MODE
  PARAMETERS$SAVE_INTERMEDIATE=SAVE_INTERMEDIATE
  PARAMETERS$TXDB=TXDB
  
  # first source the c++ fiel
  # sourceCpp(system.file('extdata','Comp.cpp',package='MeTPeak'))
  # check annotation
  if (is.na(PARAMETERS$GENOME) & is.na(PARAMETERS$GENE_ANNO_GTF)) { 
    stop("must specify the genome assembly or provide a gene gtf file for MeTPeak to work!", 
         call. = TRUE, domain = NULL)}
  
  # dependent variables
  if (is.na(PARAMETERS$READ_LENGTH)) {PARAMETERS$READ_LENGTH=.get.bam.read.length(PARAMETERS$IP_BAM[1])}
  if (is.na(PARAMETERS$MINIMAL_PEAK_LENGTH)) {PARAMETERS$MINIMAL_PEAK_LENGTH=PARAMETERS$FRAGMENT_LENGTH/2}
  if (is.na(PARAMETERS$PEAK_CUTOFF_PVALUE)) {PARAMETERS$PEAK_CUTOFF_TYPE="FDR"} else  {PARAMETERS$PEAK_CUTOFF_TYPE="PVALUE"}
  if (is.na(PARAMETERS$OUTPUT_DIR)) {PARAMETERS$OUTPUT_DIR=getwd()}
  
  # algrithm ##################################################
  
  # read gene annotation
  ANNOTATION = .read.gtf(PARAMETERS)
  ANNOTATION_BATCH_ID = .divide.anno.into.batches(ANNOTATION)
  
  # index bam files
  # get bam file
  bam=c(PARAMETERS$IP_BAM,PARAMETERS$INPUT_BAM)
  no_bam_files=length(bam)
  for (ibam in 1:no_bam_files) {file=bam[ibam]
  if (! file.exists(paste(file,'.bai',sep=""))) {
    print(paste("Stage: index bam file", file))
    indexBam(file)
  }}
  
  # get reads count report
  SAMPLE_ID = .get.sample.id(PARAMETERS)
  BAM_CHRS = .get.bam.chrs(PARAMETERS$IP_BAM[1])
  
  # get reads count
  print("Get Reads Count ...")
  print("This step may take a few hours ...")
  if (is.na(PARAMETERS$TESTING_MODE)) {no_batch_to_run=max(ANNOTATION_BATCH_ID)} else {no_batch_to_run=PARAMETERS$TESTING_MODE}
  noGroups=ceiling(no_batch_to_run/170);
  group_bar=round(seq(from = 1, to = no_batch_to_run +1,length.out=noGroups+1))
  # get space
  #READS_COUNT=.get.reads.count(1,PARAMETERS,ANNOTATION,ANNOTATION_BATCH_ID,BAM_CHRS,no_bam_files,bam)
  #READS_COUNT=READS_COUNT[integer(0),]
  #READS_COUNT_FORMAT=READS_COUNT
  READS_COUNT=list()
  # groups
  for (igroup in 1:noGroups){
    print( paste(as.character(signif(igroup*100/noGroups, digits = 3)),"%") )
    temp=list()
    batch_batch=group_bar[igroup]:(group_bar[igroup+1]-1)
    #reads_count_group=READS_COUNT_FORMAT
    no_batch=length(batch_batch)
    for (ibatch in 1:no_batch){
      temp[[ibatch]]=my.get.reads.count(batch_batch[ibatch],PARAMETERS,ANNOTATION,ANNOTATION_BATCH_ID,BAM_CHRS,no_bam_files,bam)
      #temp[[ibatch]]=.get.reads.count(batch_batch[ibatch],PARAMETERS,ANNOTATION,ANNOTATION_BATCH_ID,BAM_CHRS,no_bam_files,bam)
    }
    #for (ibatch in 1:no_batch){
    #  reads_count_group=rbind(reads_count_group,temp[[ibatch]])
    #}
    READS_COUNT[[igroup]]=dplyr::bind_rows(temp)
    #READS_COUNT=rbind(READS_COUNT,reads_count_group)
  }
  READS_COUNT = dplyr::bind_rows(READS_COUNT)
  saveRDS(object = READS_COUNT,file = paste0(PARAMETERS$OUTPUT_DIR,"READS_COUNT.rds"))
}