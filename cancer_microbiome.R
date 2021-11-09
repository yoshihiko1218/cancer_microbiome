FilePrepare <- function(CancerMiExp,NormalMiExp,GeneExp,clinical,vars,out="./",IAGenes="./IAGene.txt",fraction=0.1,ifFilter='or',acomp=TRUE,adjplot='adj',genePick = 30) {
  #Prepare Packages
  packages <- c("tidyverse","compositions","varhandle","plyr","progress","taxize")
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))  
  }
  #biocPackage <- c("DESeq2")
  #if (length(setdiff(biocPackage, rownames(installed.biocPackage()))) > 0) {
  #  BiocManager::install(setdiff(biocPackage, rownames(installed.biocPackage())))  
  #}
  library(DESeq2)
  library(tidyverse)
  library(compositions)
  library(plyr)
  library(DESeq2)
  library(progress)
  library(taxize)
  Sys.setenv(ENTREZ_KEY = "42153280b1c121b2fe7098ea2111a8c30d08")
  
  YoshiiPack <- list()
  
  cli <- read.csv(clinical,header=FALSE,row=1,sep='\t')
  YoshiiPack$clinical <- cli
  varlist <- read.csv(vars,header=TRUE,sep='\t')
  YoshiiPack$vars <- varlist
  #Read IA gene list
  IAGene <- read.csv(IAGenes,header=FALSE)
  
  #Reading date of sequence
  sub <- c()
  dateC<-read.csv(CancerMiExp,header=FALSE, row.names = 1,sep='\t',nrows = 2)
  for(i in colnames(dateC)){
    sub <- c(sub,substr(dateC[i][[1]][1],10,15))
  }
  dateC[1,] <- sub
  dateC <- t(dateC)
  colnames(dateC)[1] <- "Date"
  
  sub <- c()
  dateN<-read.csv(NormalMiExp,header=FALSE, row.names = 1,sep='\t',nrows = 2)
  for(i in colnames(dateN)){
    sub <- c(sub,substr(dateN[i][[1]][1],10,15))
  }
  dateN[1,] <- sub
  dateN <- t(dateN)
  colnames(dateN)[1] <- "Date"
  
  dateComb <- rbind(dateC, dateN)
  YoshiiPack$date <- dateComb
  
  #Reading and Prepare Cancer Mic Expression
  MicC<-read.csv(CancerMiExp,header=FALSE,skip = 1, row.names = 1,sep='\t')
  MicCTrans <- t(MicC)
  MicCTrans <- as_tibble(MicCTrans)
  names(MicCTrans) <- gsub(x = names(MicCTrans), pattern = "\\|", replacement = ".")
  MicCTrans <- data.frame(append(MicCTrans, c(Type='Cancer'), after=1))

  #Reading and Prepare Normal Mic Expression
  MicN<-read.csv(NormalMiExp,header=FALSE,skip = 1,row.names = 1,sep='\t')
  MicNTrans <- t(MicN)
  MicNTrans <- as_tibble(MicNTrans)
  names(MicNTrans) <- gsub(x = names(MicNTrans), pattern = "\\|", replacement = ".")
  MicNTrans <- data.frame(append(MicNTrans, c(Type='Normal'), after=1))

  #Reading and Prepare Gene Expression
  Gene<-read.csv(GeneExp,header=FALSE,sep='\t')
  design <- t(Gene[1,][2:length(Gene[1,])])
  colnames(design) <- "Type"
  colnames(Gene) <- Gene[2,]
  Gene <- Gene[-1,]
  Gene <- Gene[-1,]
  rownames(Gene) = make.names(Gene$Patient, unique=TRUE)
  Gene <- Gene[,-1]
  Gene <- replace(Gene, Gene == "", 0)
  dds <- DESeqDataSetFromMatrix(data.matrix(Gene), DataFrame(design),~Type)
  dds <- DESeq(dds)
  res <- results(dds)
  write.csv( as.data.frame(res), file=paste0(out,"DEresults.txt"))
  signames <- rownames(subset(res[order(res$padj),], padj < 0.05))
  if(length(intersect(rownames(subset(res[order(res$padj),], padj < 0.05)),IAGene["V1"][,1])) >= genePick){
    criticalGene <- intersect(rownames(subset(res[order(res$padj),], padj < 0.05)),IAGene["V1"][,1])[1:genePick]
  }
  else{
    criticalGene <- intersect(rownames(subset(res[order(res$padj),], padj < 0.05)),IAGene["V1"][,1])
  }
  
  bothGene<-read.csv(GeneExp,header=TRUE,sep='\t')
  Gene <- cbind(bothGene['X'],bothGene[, grepl("^Cancer", names(bothGene))])
  rownames(Gene) = make.names(Gene$X, unique=TRUE)
  Gene<-Gene[,-1]
  YoshiiPack$GeneExp <- Gene
  #Gene<-read.csv(GeneExp,header=FALSE,row.names = 1,sep='\t')
  GeneTrans <- t(Gene)
  GeneTrans <- as_tibble(GeneTrans)
  
  
  comb <- rbind.fill(MicCTrans,MicNTrans)
  comb[is.na(comb)] <- 0
  comb <- comb
  
  #Remove By Fraction for Cancer Mic
  sum <- colSums(comb != 0)
  begin <- comb[,1:2]
  for(i in 3:ncol(comb)){
    if(sum[i]>(fraction*nrow(comb))){
      begin <- cbind(begin,comb[i])
    }
  }
  comb <- begin
  
  
  
  if(acomp == TRUE){
    comb<-cbind(comb[,1:2],as_tibble(acomp(comb[,3:ncol(comb)])))
  }
  else{
    comb<-cbind(comb[,1:2],as_tibble(comb[,3:ncol(comb)]))
  }
  #message(comb[,1])
  comb <- as_tibble(comb)
  BiMic <- Bi(comb[1:length(row.names(MicCTrans)),])
  #message(BiMic[,1])
  names(BiMic) <- gsub(x = names(BiMic), pattern = "ti\\.", replacement = "")
  names(comb) <- gsub(x = names(comb), pattern = "ti\\.", replacement = "")
  newName <- colnames(BiMic)[1:2]
  oldID <- names(BiMic)
  
  taxize_oil_class <- classification(oldID[3:length(BiMic)], db = "ncbi")
  for ( i in 1:(length(colnames(BiMic))-2)){
    newName <- c(newName,taxize_oil_class[[i]][[1]][length(taxize_oil_class[[i]][[1]])])
  }
  newName <- make.unique(newName, sep=" ")
  newName <- gsub(x = newName, pattern = " ", replacement = "")
  newName <- gsub(x = newName, pattern = "\\.", replacement = "")
  newName <- gsub(x = newName, pattern = "-", replacement = "")
  newName <- gsub(x = newName, pattern = "\\[", replacement = "")
  newName <- gsub(x = newName, pattern = "\\]", replacement = "")
  newName <- gsub(x = newName, pattern = "\\+", replacement = "")
  newName <- gsub(x = newName, pattern = "\\'", replacement = "")
  newName <- gsub(x = newName, pattern = "\\/", replacement = "")
  newName <- gsub(x = newName, pattern = "\\#", replacement = "")
  names(BiMic) <- newName
  names(comb) <- newName
  
  YoshiiPack$BiMic <- BiMic
  YoshiiPack$comb <- comb
  
  common<- intersect(names(GeneTrans),criticalGene)
  GeneTrans <- cbind(GeneTrans[,1],GeneTrans[,common])
  
  if(ifFilter == "Surv" || ifFilter == "Diff" || ifFilter == 'and' || ifFilter == 'or'){
    
    if(ifFilter == 'and' || ifFilter == "Diff" || ifFilter == 'or'){
      #run Differental Abundance
      dir.create(file.path("./DifferentalAbundance"), showWarnings = FALSE)
      goodMicDE <- Kruskal(task="DifferentialAbundance",out="./DifferentalAbundance/",matinput = cbind(cbind(comb[1],comb[,3:length(names(comb))]),comb["Type"]), CatNum = 1,NumNum = length(names(comb))-2,ifplot=TRUE,cutoff=0.1,adjplot=adjplot)

    }
    
    if(ifFilter == 'and' || ifFilter == "Surv" || ifFilter == 'or'){
      #runSurvival
      dir.create(file.path("./Survival"), showWarnings = FALSE)
      goodMicSurv <- Survival(df =YoshiiPack,ifReturn = T,clinical,out="./Survival/")

    }
    
    if(ifFilter == "Surv"){
      if(length(goodMicSurv)==0){
        message("Filter does not work for this, please run without filter(By Survival)")
      }
      else{
        pickup <- goodMicSurv
        YoshiiPack$comb <- cbind(comb[,1:2],comb[pickup])
        YoshiiPack$BiMic <- cbind(BiMic[,1:2],BiMic[pickup])
      }
    }
    else if(ifFilter == "Diff"){
      if(length(goodMicDE)==0){
        message("Filter does not work for this, please run without filter(By DiffAbundance)")
      }
      else{
        pickup <- goodMicDE
        YoshiiPack$comb <- cbind(comb[,1:2],comb[pickup])
        YoshiiPack$BiMic <- cbind(BiMic[,1:2],BiMic[pickup])
      }
    }
    else if(ifFilter == 'and'){
      pickup <- intersect(goodMicSurv,goodMicDE)
      if(length(pickup)==0){
        message("Filter does not work because there's no intersection!! : (")
      }
      else{
        YoshiiPack$comb <- cbind(comb[,1:2],comb[pickup])
        YoshiiPack$BiMic <- cbind(BiMic[,1:2],BiMic[pickup])
      }
    }
    else if(ifFilter == 'or'){
      pickup <- unique(c(goodMicSurv,goodMicDE))
      if(length(pickup)==0){
        message("Filter does not work because there's no significant Microbe!?!? : (")
      }
      else{
        YoshiiPack$comb <- cbind(comb[,1:2],comb[pickup])
        YoshiiPack$BiMic <- cbind(BiMic[,1:2],BiMic[pickup])
      }
    }
  }
  else{
    YoshiiPack$BiMic <- BiMic
    YoshiiPack$comb <- comb
  }
  GeneMic <- merge(GeneTrans, YoshiiPack$BiMic, by="Patients", all=F,suffixes="Patients")
  GeneMic <- GeneMic %>% filter(Patients != "Patients")
  YoshiiPack$GeneMic <- subset(GeneMic, select=-c(Type))
  YoshiiPack$GeneTrans <- GeneTrans
  return(YoshiiPack)
}

DifferentalAbundance <- function(df,out = "./DifferentalAbundance/", cutoff = 0.1,adjplot='adj'){
  dir.create(file.path(out), showWarnings = FALSE)
  Kruskal(task="IndividualDifferentialAbundance",out=out,matinput = cbind(cbind(df$comb[1],df$comb[,3:length(names(df$comb))]),df$comb["Type"]), CatNum = 1,NumNum = length(names(df$comb))-2,ifplot=TRUE,cutoff=0.1,adjplot=adjplot)
}

IAKruskal <- function(df,out="./IACorrelation/",cutoff=0.1,adjplot = 'adj'){
  dir.create(file.path(out), showWarnings = FALSE)
  Kruskal(task="IA",out=out,matinput=df$GeneMic,CatNum = length(names(df$BiMic))-2, NumNum = length(names(df$GeneTrans))-1,cutoff=cutoff,adjplot=adjplot)
}

Bi <- function(CorMic){
  BiMic <- CorMic[,1:2]
  MicMatrix <- as.matrix(CorMic)
  for(i in 3:ncol(CorMic)){
    colMean <- mean(as.numeric(MicMatrix[,i]))
    newCol <- c()
    for(j in 1:nrow(CorMic)){
      if(as.numeric(MicMatrix[,i][j])>=colMean){
        newCol <- c(newCol,"HIGH")
      }
      else{
        newCol <- c(newCol,"LOW")
      }
    }
    #message(paste0('mean:',newCol))
    BiMic[names(CorMic)[i]]<-newCol
  }
  return(BiMic)
}


Kruskal <- function(task,out,matinput,CatNum, NumNum,ifplot = TRUE, cutoff=0.1,adjplot){
  message(">>> Running Kruskal")
  pb <- progress_bar$new(
    format = "  Kruskal [:bar] :percent eta: :eta",
    total = CatNum*NumNum, clear = FALSE, width= 60)
  goodPair <- c()
  mat <- matrix(, nrow = CatNum*NumNum+1, ncol = 5)
  mat[1,1] <- "Continuous Variable"
  mat[1,2] <- "Categorical Variable"
  mat[1,3] <- "Pvalue"
  mat[1,4] <- "ifSignicicant"
  mat[1,5] <- "log2FoldChange"
  n = 1
  ColNames <- names(matinput)
  CoMatrix <<- as.matrix(matinput)
  goodMic<-c()
  pb$tick(0)
  for (i in 2:(NumNum+1)) {
    for (j in (NumNum+2):(NumNum+1+CatNum)){
      pb$tick()
      Sys.sleep(1 / 100)
      NonNAindex <- which(!is.na(CoMatrix[,j]))
      if(length(NonNAindex) == 0){
        firstNonNA = 1
      }
      else{
        firstNonNA <- min(NonNAindex)
      }
      if(!(all(CoMatrix[,j]!=CoMatrix[,j][firstNonNA]))){
        tryCatch({
        krs <- kruskal.test(CoMatrix[,i]~CoMatrix[,j])
        n = n+1
        mat[n,1] <- ColNames[i]
        mat[n,2] <- ColNames[j]
        mat[n,3] <- krs$p.value
        if (krs$p.value <= 0.05) {
          mat[n,4] <- 'TRUE'
          goodPair <- c(goodPair,ColNames[i],ColNames[j])
          if(task=="DifferentialAbundance"){
            goodMic <- c(goodMic,ColNames[i])
          }
        } else {
          mat[n,4] <- 'FALSE'
        }
        if(task == "DifferentialAbundance" || task == "IndividualDifferentialAbundance"){
          sub1 <- dplyr::filter(as_tibble(CoMatrix), Type == "Cancer")
          sub2 <- dplyr::filter(as_tibble(CoMatrix), Type == "Normal")
          mat[n,5] <- log2(mean(as.numeric(sub1[[i]]))/(mean(as.numeric(sub2[[i]]))))
        }
        else if(task == "CIBERSORT" || task == "IA"){
          #message( ColNames[j])
          sub1 <<- dplyr::filter(as_tibble(CoMatrix), ColNames[j] == "HIGH")
          sub2 <<- dplyr::filter(as_tibble(CoMatrix), ColNames[j] == "LOW")
          mat[n,5] <- log2(mean(as.numeric(sub1[[i]]))/(mean(as.numeric(sub2[[i]]))))
        }
        })
      }
    }
  }
  goodPair <- goodPair
  write.table(CoMatrix,paste0(paste0(out,task),"KWInput.txt"), sep="\t",  col.names=TRUE,row.names = FALSE,quote = FALSE)
  write.table(mat,paste0(paste0(out,task),"KWResult.txt"), sep="\t",  col.names=FALSE,row.names = FALSE,quote = FALSE)
  if(ifplot == TRUE){
    if(length(goodPair)!=0){
      Plotting(out,paste0(paste0(out,task),"KWInput.txt"),paste0(paste0(out,task),"KWResult.txt"),cutoff=cutoff,adjplot=adjplot)
    }
  }
  if(task=="DifferentialAbundance"){
    return(goodMic)
  }
}

Survival <- function(df,ifReturn=F,clinical=F,out="./Survival/"){
  packages <- c("survival")
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))  
  }
  library(survival)
  message(">>> Running Survival")
  dir.create(file.path(out), showWarnings = FALSE)
  if(clinical != F){
    cli <- read.csv(clinical,header=FALSE,row=1,sep='\t')
  }
  else{
    cli <- df$clinical
  }
  status <- cli["patient.bcr_patient_uuid",]
  status <- rbind(status,cli["patient.vital_status",])
  status <- rbind(status,cli["patient.days_to_death",])
  status <- rbind(status,cli["patient.days_to_last_followup",])
  time <- c()
  St <- c()
  for(i in 1:length(status["patient.bcr_patient_uuid",])){
    if(status["patient.vital_status",][i]=="dead"){
      St <- c(St,1)
      time <- c(time,as.numeric(status["patient.days_to_death",][i]))
    }
    else{
      St <- c(St,0)
      time <- c(time,as.numeric(status["patient.days_to_last_followup",][i]))
    }
  }
  status["status",] <- St
  status["time",] <- time
  status <- status[c("patient.bcr_patient_uuid","status","time"),]
  rownames(status)[1]<-"Patients"
  statusTrans <- as_tibble(t(status), rownames = "row_names")
  statusTrans<-statusTrans[,2:ncol(statusTrans)]
  StatusComb <- merge(statusTrans, df$BiMic, by="Patients", all=F,suffixes="Patients")
  StatusComb <- StatusComb %>% filter(Patients != "Patients")
  micnames <- names(StatusComb)
  #StatusComb <<- StatusComb
  goodMic <- c()
  options(warn=-1)
  pb <- progress_bar$new(
    format = "  Survival [:bar] :percent eta: :eta",
    total = length(micnames)-4, clear = FALSE, width= 60)
  pb$tick(0)
  report <- "MicrobeName\tCorrelaton\tPval\tifSignificant\n"
  for(i in 5:length(micnames)){
    pb$tick()
    Sys.sleep(1 / 100)
    #message(micnames[i])
    iformula <- as.formula(sprintf("Surv(as.numeric(time), as.numeric(status)) ~ %s", micnames[i]))
    savecox <-summary(coxph(iformula, data=StatusComb))
    pVal <- savecox[[7]][1,5]
    if(savecox[[7]][2] > 1){
      corre <- "negative"
    }
    else{
      corre <- "positive"
    }
    report <- paste0(paste0(paste0(paste0(paste0(paste0(report,micnames[i]),'\t'),corre),'\t'),pVal),'\t')
    if(pVal <= 0.05){
      report <- paste0(paste0(report,"TRUE"),'\n')
    }
    else{
      report <- paste0(paste0(report,"FALSE"),'\n')
    }
    if(pVal <= 0.05){
      goodMic <- c(goodMic,micnames[i])
      pdf(paste0(out,paste0(micnames[i],"Survival.pdf")))
      plot(survfit(iformula,data=StatusComb), main=micnames[i], col=c("red","blue"), xlab="Time (Days)", ylab="Proportion Surviving", mark.time=TRUE)
      legend("topright",1.0,c("High Expression","Low Expression"), lty=c(1,1), col=c("red","blue"))
      dev.off()
    }
  }
  options(warn=0)
  fileConn<-file(paste0(out,"Report.txt"))
  writeLines(report, fileConn)
  close(fileConn)
  if(ifReturn==T){
    return(goodMic)
  }
}

ciber <- function(df, Gene=F,Rfile="./CIBERSORT.R",Target="./LM22.txt",cutoff=0.1,out="./CIBERSORT/",adjplot='adj'){
  dir.create(file.path(out), showWarnings = FALSE)
  packages <- c("remotes")
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))  
  }
  #remotes::install_github("icbi-lab/immunedeconv")
  library(immunedeconv)
  #GeneExp<-read.csv(Gene,header=FALSE,row.names = 1,sep='\t')
  if(Gene != F){
  bothGene<-read.csv(Gene,header=TRUE,sep='\t')
  GeneExp <- cbind(bothGene['X'],bothGene[, grepl("^Cancer", names(bothGene))])
  rownames(GeneExp) = make.names(GeneExp$X, unique=TRUE)
  GeneExp<-GeneExp[,-1]
  }
  else{
    GeneExp <- YoshiiPack$GeneExp
  }
  
  colnames(GeneExp)<-GeneExp[1,]
  GeneExp <- GeneExp[2:length(row.names(GeneExp)),]
  GeneExp[is.na(GeneExp)] <- 0
  GeneExp <- replace(GeneExp, GeneExp == "", 0)
  set_cibersort_mat(Target)
  set_cibersort_binary(Rfile)
  ImmuneCells <- deconvolute(GeneExp,"cibersort", tumor = TRUE)
  message(">>> cibersort Done Running")
  row.names(ImmuneCells) <- as.character(unlist(ImmuneCells[,1]))
  ImmuneTrans <- t(ImmuneCells)
  ImmuneTrans <- ImmuneTrans[-1,]
  ImmuneTrans <- as.data.frame(ImmuneTrans)
  ImmuneTrans <- rownames_to_column(ImmuneTrans,"Patients")
  CombinedImmune <- merge(ImmuneTrans, df$BiMic, by="Patients", all=F,suffixes="Patients")
  CombinedImmune <- subset(CombinedImmune, select=-c(Type))
  Kruskal(task="CIBERSORT",out = out,matinput = CombinedImmune,CatNum = (ncol(df$BiMic)-2),NumNum = (ncol(ImmuneTrans)-1),cutoff= cutoff,adjplot=adjplot)
}

Plotting <- function(out,KWInput,KWResult,cutoff,adjplot){
  message(">>> Plotting")
  inputmat <- read.csv(KWInput, header=FALSE,sep='\t')
  colnames(inputmat) <- inputmat[1,]
  inputmat <- inputmat[-1,]
  inputmat <<- as_tibble(inputmat)
  res <- read.csv(KWResult, header=FALSE,sep='\t')
  res <- as_tibble(res)
  res <- res %>% filter(V4 == "TRUE")
  goodPair <- c()
  for(i in 1:length(row.names(res))){
    goodPair <- c(goodPair,res["V1"][i,])
    goodPair <- c(goodPair,res["V2"][i,])
  }
  #goodPair<-goodPair
  pb <- progress_bar$new(
    format = "  Plotting [:bar] :percent eta: :eta",
    total = length(goodPair)/2, clear = FALSE, width= 60)
  for(i in 1:(length(goodPair)/2)){
    pb$tick()
    Sys.sleep(1 / 100)
    gene = goodPair[[2*i-1]]
    variable = goodPair[[2*i]]
    plotStats<-boxplot(as.numeric(inputmat[,gene])~as.character(inputmat[,variable]),plot=FALSE,na.rm = TRUE)
    #find number of categories
    ncat<-ncol(plotStats$stats)
    
    #find the median of each box
    mat<-matrix(,nrow=ncat,ncol=1)
    for(i in 1:ncat){ 
      mat[i,]<-plotStats$stats[3,i] #always five rows in plotStats, 3rd row is median
    }
    
    #find differences between medians
    diffMat<-diff(mat)
    
    #check for trend 
    trendCheck<-0 #initialize trend variable
    for (i in 1:(ncat-1)){
      trendCheck<-trendCheck+sign(diffMat[i,]) 
    } #add signs of differences to determine if trend exists
    #message(trendCheck)
    # if there is trend, check signifcance of median differences
    if (abs(trendCheck)>=1){ # give allowance of 1 category in opposite direction
      weight<-matrix(,nrow=(ncat-1),ncol=1)
      for (i in 1:(ncat-1)) {
        weight[i,]<-diffMat[i,]/max(as.numeric(inputmat[,gene]),na.rm = TRUE) #find difference between medians as percentage of max value
      }
      #message(weight)
      if(!is.na(weight)){
        if(abs(mean(weight,na.rm=TRUE))> cutoff){
          if(adjplot == 'nonadj' || adjplot == 'both'){
            pdf(paste0(out,paste0(gene, variable,".pdf", sep="")))
            boxplot(as.numeric(inputmat[,gene])~inputmat[,variable], xlab=variable ,ylab=gene, outline=TRUE,las=2,na.rm=TRUE)
            dev.off()
          }
          if(adjplot == 'both' || adjplot == 'adj'){
            outliers <- boxplot(as.numeric(inputmat[,gene]), plot=FALSE)$out
            if(is_empty(outliers) == FALSE ){
              mid<-inputmat
              mid<- mid[-which(as.numeric(mid[,gene]) %in% outliers),]
              pdf(paste0(out,paste0(gene, variable,"adjusted.pdf", sep="")))
              boxplot(as.numeric(mid[,gene])~mid[,variable], xlab=variable ,ylab=gene, outline=TRUE,las=2,na.rm=TRUE)
              dev.off()
            }
          }
        }
      }
    } 
  }
}


ClinicalVar <- function(df,out="./ClinicalVar/",vars = F,clinical=F,cutoff=0.1,adjplot="adj"){
  packages <- c("varhandle")
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))  
  }
  library(varhandle)
  if(clinical != F){
    cli <- read.csv(clinical,header=FALSE,row=1,sep='\t')
  }
  else{
    cli <- df$clinical
  }
  if(vars != F){
    varlist <- read.csv(vars,header=TRUE,sep='\t')
  }
  else{
    varlist <- df$vars
  }
  dir.create(file.path(out), showWarnings = FALSE)
  headNum <- cli[names(varlist)[1],]
  headCat <- cli[names(varlist)[1],]
  used <- intersect(row.names(cli),varlist[,1])
  for(i in used){
    NonNAindex <- which(!is.na(cli[i,]))
    if(length(NonNAindex) == 0){
      firstNonNA = 1
    }
    else{
      firstNonNA <- min(NonNAindex)
    }
    if(firstNonNA != Inf){
      if(check.numeric(cli[i,][[firstNonNA]])==TRUE){
        headNum <- rbind(headNum,cli[i,])
      }
      else{
        headCat <- rbind(headCat,cli[i,])
      }
    }
  }
  saverow <- rownames(headCat)
  headCat <- data.frame(lapply(headCat, function(x) {gsub("mx", NA, x)}))
  headCat <- data.frame(lapply(headCat, function(x) {gsub("gx", NA, x)}))
  headCat <- data.frame(lapply(headCat, function(x) {gsub("nx", NA, x)}))
  headCat <- data.frame(lapply(headCat, function(x) {gsub("tx", NA, x)}))
  headCat <- data.frame(lapply(headCat, function(x) {gsub("t2a", "t2", x)}))
  headCat <- data.frame(lapply(headCat, function(x) {gsub("t2b", "t2", x)}))
  headCat <- data.frame(lapply(headCat, function(x) {gsub("t3a", "t3", x)}))
  headCat <- data.frame(lapply(headCat, function(x) {gsub("t3b", "t3", x)}))
  headCat <- data.frame(lapply(headCat, function(x) {gsub("t4a", "t4", x)}))
  headCat <- data.frame(lapply(headCat, function(x) {gsub("t4b", "t4", x)}))
  
  #headCat <- data.frame(lapply(headCat, function(x) {gsub("NA", NA, x)}))
  #headCat <<- headCat
  #headCat <- data.frame(lapply(headCat, function(x) {gsub(NA, '', x)}))
  
  rownames(headCat) <- saverow
  rownames(headNum)[1]<-"Patients"
  rownames(headCat)[1]<-"Patients"
  NumTrans <- as_tibble(t(headNum), rownames = "row_names")
  NumTrans <- NumTrans[,2:ncol(NumTrans)]
  CatTrans <- as_tibble(t(headCat), rownames = "row_names")
  CatTrans <-CatTrans[,2:ncol(CatTrans)]
  CatTrans <- CatTrans
  NumCombined <- merge(NumTrans, df$BiMic, by="Patients", all=F,suffixes="Patients")
  NumCombined <- subset(NumCombined, select=-c(Type))
  NumCombined <- NumCombined %>% filter(Patients != "Patients")
  NumCombined <- NumCombined
  Kruskal(task = "Clinical", out = paste0(out,"num"),matinput = NumCombined,CatNum = (ncol(df$BiMic)-2),NumNum = (ncol(NumTrans)-1),cutoff = cutoff,adjplot=adjplot)
  CatCombined <- merge(df$comb[1:length(row.names(df$BiMic)),],CatTrans, by="Patients", all=F,suffixes="Patients")
  CatCombined <- subset(CatCombined, select=-c(Type))
  CatCombined <- CatCombined %>% filter(Patients != "Patients")
  CatCombined <- CatCombined
  Kruskal(task = "Clinical",out = paste0(out,"cat"),matinput = CatCombined,CatNum = ncol(CatTrans)-1,NumNum = ncol(df$BiMic)-2,cutoff = cutoff,adjplot=adjplot)
}

RunGSEA <- function(GeneExp,skipInput=F,MicrobeList = F,prefix="test",gmtfile ="./c2c6c7_combined.gmt",input.chip = "NOCHIP", gene.ann = "",  
                      gs.ann = "", output.directory = paste0(getwd(),'/'), reshuffling.type = "sample.labels", 
                      nperm = 1000, weighted.score.type = 1, nom.p.val.threshold = -1, fwer.p.val.threshold = -1, 
                      fdr.q.val.threshold = 0.25, topgs = 20, adjust.FDR.q.val = F, gs.size.threshold.min = 15, 
                      gs.size.threshold.max = 500, reverse.sign = F, preproc.type = 0, random.seed = as.integer(Sys.time()), 
                      perm.type = 0, fraction = 1, replace = F, collapse.dataset = FALSE, collapse.mode = "NOCOLLAPSE", 
                      save.intermediate.results = F, use.fast.enrichment.routine = T, gsea.type = "GSEA", 
                      rank.metric = "S2N"){
  packages <- c("plyr","varhandle","compositions","tidyverse")
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))  
  }
  library(plyr)
  library(varhandle)
  library(compositions)
  library(tidyverse)
  dir.create(file.path(output.directory), showWarnings = FALSE)
  #message(paste0(output.directory,"Inputs"))
  dir.create(file.path(paste0(output.directory,"Inputs")), showWarnings = FALSE)
  dir.create(file.path(paste0(output.directory,"Results")), showWarnings = FALSE)
  
  dir.create(file.path(output.directory), showWarnings = FALSE)
  #AllGene<-read.csv(GeneExp,header=FALSE,row.names = 1,sep='\t')
  bothGene<-read.csv(GeneExp,header=TRUE,sep='\t')
  AllGene = cbind(bothGene['X'],bothGene[, grepl("^Cancer", names(bothGene))])
  rownames(AllGene) = make.names(AllGene$X, unique=TRUE)
  AllGene<-AllGene[,-1]
  AllGene <- replace(AllGene, AllGene == "", 0)
  AllGeneTrans <- t(AllGene)
  AllGeneTrans <- as_tibble(AllGeneTrans)
  
  if(skipInput == F){
    makeGCT(output.directory,prefix,AllGeneTrans)
    makeCLS(prefix,AllGeneTrans,gmt,output.directory)
  }
  if(MicrobeList == F){
    MicList <- names(BiMic)[3:length(BiMic[1,])]
  }
  else{
    MicList <- read.csv(MicrobeList,header=F,sep="\t")
    MicList <- MicList["V1"][,1]
  }
  MicList <- names(BiMic)[3:length(BiMic[1,])]
  HereGSEA(paste0(paste0(output.directory,"Inputs/"),paste0(prefix,".gct")), MicList,input.chip = input.chip, gene.ann = gene.ann, gmt = gmtfile, 
          gs.ann = gs.ann, out = output.directory, reshuffling.type = reshuffling.type, 
          nperm = nperm, weighted.score.type = weighted.score.type, nom.p.val.threshold = nom.p.val.threshold, fwer.p.val.threshold = fwer.p.val.threshold, 
          fdr.q.val.threshold = fdr.q.val.threshold, topgs = topgs, adjust.FDR.q.val = adjust.FDR.q.val, gs.size.threshold.min = gs.size.threshold.min, 
          gs.size.threshold.max = gs.size.threshold.max, reverse.sign = reverse.sign, preproc.type = preproc.type, random.seed = random.seed, 
          perm.type = perm.type, fraction = fraction, replace = replace, collapse.dataset = collapse.dataset, collapse.mode = collapse.mode, 
          save.intermediate.results = save.intermediate.results, use.fast.enrichment.routine = use.fast.enrichment.routine, gsea.type = gsea.type, 
          rank.metric = rank.metric)
}

HereGSEA <- function(gctName,clsList,input.chip, gene.ann, gmt, 
                    gs.ann, out, reshuffling.type, 
                    nperm, weighted.score.type, nom.p.val.threshold, fwer.p.val.threshold, 
                    fdr.q.val.threshold, topgs, adjust.FDR.q.val, gs.size.threshold.min, 
                    gs.size.threshold.max, reverse.sign, preproc.type, random.seed, 
                    perm.type, fraction, replace, collapse.dataset, collapse.mode, 
                    save.intermediate.results, use.fast.enrichment.routine, gsea.type, 
                    rank.metric){
  packages <- c("devtools")
  if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))  
  }
  library(devtools)
  #install_github("GSEA-MSigDB/GSEA_R")
  library(GSEA)
  message(">>> Start running GSEA")
  ds <- read.csv(gctName,sep="\t",header = F)
  ds <- ds[-c(1:2), ]
  colnames(ds) <- ds[c(1), ]
  ds <- ds[-c(1), ]
  clsList <- c("ConiophoraputeanaRWD64598SS2")
  for(i in clsList){
    #message(paste0(paste0(paste0(output.directory,"Inputs/"),i),".cls"))
    #dir.create(file.path(paste0(paste0(output.directory,"Results/"),i)), showWarnings = FALSE)
    GSEA(ds, paste0(paste0(paste0(out,"Inputs/"),i),".cls"), gs.db = gmt, output.directory = paste0(paste0(paste0(out,"Results/"),i),"/"),doc.string = paste0(i,"gsea_result"), random.seed = 101)
    message(paste0(i,",done"))
  }
}

makeBASH <- function(time,email,prefix,LCfolder,bin){
  message(">>> creating bash file")
  bash <- matrix(, nrow = 18, ncol = 1)
  bash[1,1] <- "#!/bin/bash"
  bash[2,1] <- paste0(paste0("#SBATCH --job-name=\"",prefix),"_GSEA\"")
  bash[3,1] <- paste0(paste0("#SBATCH --output=\"",prefix),"_GSEA.out\"")
  bash[4,1] <- "#SBATCH -p compute"
  bash[5,1] <- "#SBATCH --nodes=1"
  bash[6,1] <- "#SBATCH --ntasks-per-node=24"
  bash[7,1] <- "#SBATCH --export=ALL"
  bash[8,1] <- paste0(paste0("#SBATCH -t  ",time),":00:00")
  bash[9,1] <- paste0("#SBATCH --mail-user=",email)
  bash[10,1] <- "#SBATCH --mail-type=all"
  bash[11,1] <- "time("
  bash[12,1] <- sprintf("export PATH=%s:$PATH",bin)
  bash[13,1] <- "export PATH=/share/apps/compute/parallel/bin:$PATH"
  bash[14,1] <- "export SLURM_NODEFILE=`generate_pbs_nodefile`"
  bash[15,1] <- "export PARALLEL=\"--workdir . --env PATH --env LD_LIBRARY_PATH --env LOADEDMODULES --env _LMFILES_ --env MODULE_VERSION --env MODULEPATH --env MODULEVERSION_STACK --env MODULESHOME --env OMP_DYNAMICS --env OMP_MAX_ACTIVE_LEVELS --env OMP_NESTED --env OMP_NUM_THREADS --env OMP_SCHEDULE --env OMP_STACKSIZE --env OMP_THREAD_LIMIT --env OMP_WAIT_POLICY --compress\""
  bash[16,1] <- "sort -u $SLURM_NODEFILE > uniqueNodes"
  bash[17,1] <- "parallel --jobs 8 --sshloginfile uniqueNodes --workdir $PWD < commands.txt"
  bash[18,1] <- ")"
  write.table( bash,paste0(LCfolder,paste0(prefix,".bash", sep = ""), sep = ""), sep="\t",  col.names=FALSE,row.names = FALSE,na="",quote = FALSE)
}

makeCLS <- function(prefix,AllGeneTrans,gmt,output.directory){
  message(">>> creating cls file")
  Co <- merge(AllGeneTrans, BiMic, by="Patients", all=F,suffixes="Patients")
  Co <- subset((Co), select=-c(Type))
  clshead <- matrix(c(length(Co[,1]),"#",2,"LOW",1,"HIGH"),nrow=2,ncol=3)
  clshead <- as.data.frame(clshead)
  micnames <- names(BiMic)
  TransCo <- t(Co)
  for(i in 3:length(BiMic)){
    nm <- matrix(micnames[i],nrow=1,ncol=1)
    nm <- as.data.frame(nm)
    va <- matrix(TransCo[micnames[i],],nrow=1,ncol=length(TransCo[micnames[i],]))
    va <- as.data.frame(va)
    clshead <- rbind.fill(clshead,rbind.fill(nm,va))
  }
  makeCLSFiles(clshead,prefix,gmt,output.directory)
  write.table( clshead,paste0(paste0(output.directory,"Inputs/"),paste0(prefix,".cls", sep = ""), sep = ""), sep="\t",  col.names=FALSE,row.names = FALSE,na="")
}


makeCLSFiles <- function(cls,prefix,gmt,output.directory){
  message(">>> creating individual cls files")
  OG <- cls
  ROWNUM<-nrow(OG)
  COLNUM<-ncol(OG)
  nmicrobes=(ROWNUM -2)*.5
  commands <- matrix(, nrow = nmicrobes, ncol = 1)
  
  x=3
  y=4
  for (i in 1:nmicrobes) {
    #message(commands)
    indiv<-matrix(,nrow=3,ncol=COLNUM)
    indiv[1,1]<-toString(OG[1,1])
    indiv[2,1]<-toString(OG[2,1])
    indiv[1,2]<-toString(OG[1,2])
    indiv[2,2]<-toString(OG[2,2])
    indiv[1,3]<-toString(OG[1,3])
    indiv[2,3]<-toString(OG[2,3])
    for (j in 1:COLNUM){
      indiv[3,j]<-toString(OG[y,j])
    }
    filename=toString(OG[x,1])
    filename=paste0(filename,".cls" ,sep = "")
    write.table(indiv, file=paste0(paste0(output.directory,"Inputs/"),filename,sep=""), sep="\t", quote=FALSE,row.names = FALSE, col.names = FALSE, na= "",)
    x=x+2
    y=y+2
    
    #message(i)
    #commands[i,1] <- sprintf("./gsea-cli.sh GSEA -res %s%s.gct -cls %s%s -gmx %s%s -collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 -permute phenotype -rnd_type no_balance -scoring_scheme weighted -rpt_label %s_antitumor_C1 -metric Pearson -sort real -order descending -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out %s",SCfolder,prefix,SCfolder,filename,SCfolder,gmt,filename,SCfolder)
  }
  #write.table(commands, file=paste0(paste0(output.directory,"/Inputs/"),paste0(prefix,"commands.txt", sep = ""), sep = ""),row.names=F, col.names=F)
}

makeGCT <- function(output.directory,prefix,AllGeneTrans){
  message(">>> creating gct file")
  Co <- merge(AllGeneTrans, BiMic, by="Patients", all=F,suffixes="Patients")
  Co <- subset((Co), select=-c(Type))
  A <- matrix(c("#1.2"), nrow=1,  ncol=1)
  A <- as.data.frame(A)
  B <- matrix(c(length(AllGeneTrans)-1,length(row.names(Co))), nrow=1,  ncol=2) 
  B <- as.data.frame(B)
  #Co <- merge(GeneTrans, BiMic, by="Patients", all=F,suffixes="Patients")
  top <- rbind.fill(A,B)
  TransCo <- t(Co)
  r <- append(rep("na",length(AllGeneTrans)-1),"Description",after=0)
  C<-matrix(r,nrow=length(AllGeneTrans),ncol=1)
  D <- cbind(C,rbind(TransCo[1,],TransCo[2:length(AllGeneTrans),]))
  D <- as.data.frame(D)
  D <- rownames_to_column(D)
  D[1,1] <- "NAME"
  names(D)[1:3] <- c("V1","V2","Random")
  p <- c()
  for(i in 1:(length(D)-2)){
    p <- c(p,i)
  }
  D[1,3:length(D)]<-BiMic[['Patients']]
  gct <- rbind.fill(top,D)
  write.table( gct,paste0(paste0(output.directory,"Inputs/"),paste0(prefix,".gct", sep = ""), sep = ""), sep="\t",  col.names=FALSE,row.names = FALSE,na="")
}

REVEALERInput <- function(CNV = F,gisticFile, mutFile,clinical,prefix="REVEALER",LCfolder="./",SCfolder="./",time=48,email="test@ucsd.edu"){
  dir.create(file.path(LCfolder), showWarnings = FALSE)
  if(CNV == F){
  inputMaker(gisticFile, mutFile,LCfolder)
  }
  else{
    CNVdata <- read.delim(CNV,header=F)
  }
  clin <- read.delim(clinical,header =F,row.names = 1)
  matching <- clin[c("patient.samples.sample.portions.portion.analytes.analyte.aliquots.aliquot.bcr_aliquot_barcode","patient.bcr_patient_uuid"),]
  colnames(CNVdata) <- tolower(CNVdata[3,])
  colnames(matching) <- tolower(matching[1,])
  inte <- intersect(colnames(matching),colnames(CNVdata))
  newCNV <- cbind(CNVdata[,1:2],CNVdata[,inte])
  colnames(newCNV) <- c(colnames(CNVdata)[1:2],tolower(matching[,inte][2,]))
  CorMicTrans <- t(comb[1:length(row.names(BiMic)),])
  colnames(CorMicTrans) <- tolower(CorMicTrans["Patients",])
  inte2 <- intersect(colnames(CorMicTrans),colnames(newCNV))
  newCNV <<- newCNV
  finalCNV<- cbind(newCNV[,1:2],newCNV[,inte2])
  finalCNV[2,2]<-ncol(finalCNV)-2
  finalCNV[3,3:length(colnames(finalCNV))] <- colnames(finalCNV)[3:length(colnames(finalCNV))]
  write.table(finalCNV,paste0(LCfolder,paste0(prefix,"feature.gct")),sep="\t",quote=FALSE,na="",row.names=FALSE,col.names=FALSE)
  CorMicTrans <- CorMicTrans[,inte2]
  CorMicTrans <- CorMicTrans[2:length(row.names(CorMicTrans)),]
  A <- matrix(c("#1.2"), nrow=1,  ncol=1)
  A <- as.data.frame(A)
  B <- matrix(c(length(row.names(CorMicTrans))-1,length(colnames(CorMicTrans))), nrow=1,  ncol=2) 
  B <- as.data.frame(B)
  top <- rbind.fill(A,B)
  r <- append(rep("na",length(row.names(CorMicTrans))-1),"Description",after=0)
  C<-matrix(r,nrow=length(row.names(CorMicTrans)),ncol=1)
  #CorMicTrans <<- CorMicTrans
  D <- cbind(C,rbind(colnames(CorMicTrans),CorMicTrans[2:length(row.names(CorMicTrans)),]))
  #message(CorMicTrans)
  D <- as.data.frame(D)
  D <- rownames_to_column(D)
  D[1,1] <- "NAME"
  names(D)[1:3] <- c("V1","V2","Random")
  gct <- rbind.fill(top,D)
  write.table(gct, file=paste0(LCfolder,paste0(prefix,"target.gct", sep = ""), sep = ""),row.names=F, col.names=F,na="",quote=F,sep='\t')
  #bash <- matrix(, nrow = 5, ncol = 1)
  #bash[1,1] <- "#!/usr/bin/env Rscript"
  #bash[2,1] <- "args = commandArgs(trailingOnly=TRUE)"
  #bash[3,1] <- paste0(paste0("setwd(\"",SCfolder),"\")")
  #bash[4,1] <- "source(\"REVEALER_library.R\")"
  #bash[5,1] <- "switch(as.numeric(args[1]),"
  #for(i in gct[,1][4:length(gct[,1])]){
  #  sub=matrix(, nrow =12, ncol = 1)
  #  sub[1,1] <-"		{REVEALER.v1("
  #  sub[2,1] <-paste0("				  ds1                     = \"",paste0(paste0(SCfolder,prefix),"target.gct"),"\",  # Target profile dataset")
  #  sub[3,1] <-paste0("				 target.name             = \"",paste0(i,"\",       # Target profile"))
  #  sub[4,1] <-"				 target.match            = \"negative\",                    # Direction of match"
  #  sub[5,1] <-paste0("				  ds2                     = \"",paste0(paste0(SCfolder,prefix),"feature.gct"),"\",  # Features dataset")
  #  sub[6,1] <-"				 seed.names              = \"NULLSEED\",               # Seed(s) name(s)"
  #  sub[7,1] <-"				 max.n.iter              = 3,                             # Maximum number of iterations"
  #  sub[8,1] <-paste0(paste0(paste0(paste0("				 pdf.output.file          = \"Revealer_",prefix),"_"),i),".pdf\",")
  #  sub[9,1] <-"				 count.thres.low         = 3,                             # Filter out features with less than count.thres.low 1's"
  #  sub[10,1] <-"				 count.thres.high        = 409,                            # Filter out features with more than count.thres.low 1's"
  #  sub[11,1] <-"				 n.markers               = 30,                            # Number of top hits to show in heatmap for every iteration"
  #  sub[12,1] <-"			   )},"
  #  bash <- rbind(bash,sub)
  #}
  #end=matrix(, nrow =1, ncol = 1)
  #end[1][1] <- ")"
  #bash <- rbind(bash,end)
  #write.table( bash,paste0(LCfolder,paste0(prefix,".R", sep = ""), sep = ""), sep="\t",  col.names=FALSE,row.names = FALSE,na="",quote=F)
  #command = matrix(, nrow=length(gct[,1])-3,ncol=1)
  #for(i in 1:length(gct[,1])-3){
  #  command[i,1] <- paste0(paste0(LCfolder,paste0(prefix,".R ", sep = ""), sep = ""),i)
  #}
  #write.table( command,file = paste0(LCfolder,"REVEALERcommands.txt"), sep="\t",  col.names=FALSE,row.names = FALSE,na="",quote=F)
#  
  #message(">>> start making bash file for REVEALER")
  #bash <- matrix(, nrow = 17, ncol = 1)
  #bash[1,1] <- "#!/bin/bash"
  #bash[2,1] <- paste0(paste0("#SBATCH --job-name=\"",prefix),"_REVEALER\"")
  #bash[3,1] <- paste0(paste0("#SBATCH --output=\"",prefix),"_REVEALER.out\"")
  #bash[4,1] <- "#SBATCH -p compute"
  #bash[5,1] <- "#SBATCH --nodes=1"
  #bash[6,1] <- "#SBATCH --ntasks-per-node=24"
  #bash[7,1] <- "#SBATCH --export=ALL"
  #bash[8,1] <- paste0(paste0("#SBATCH -t  ",time),":00:00")
  #bash[9,1] <- paste0("#SBATCH --mail-user=",email)
  #bash[10,1] <- "#SBATCH --mail-type=all"
  #bash[11,1] <- "time("
  #bash[12,1] <- "export PATH=/share/apps/compute/parallel/bin:$PATH"
  #bash[13,1] <- "export SLURM_NODEFILE=`generate_pbs_nodefile`"
  #bash[14,1] <- "export PARALLEL=\"--workdir . --env PATH --env LD_LIBRARY_PATH --env LOADEDMODULES --env _LMFILES_ --env MODULE_VERSION --env MODULEPATH --env MODULEVERSION_STACK --env MODULESHOME --env OMP_DYNAMICS --env OMP_MAX_ACTIVE_LEVELS --env OMP_NESTED --env OMP_NUM_THREADS --env OMP_SCHEDULE --env OMP_STACKSIZE --env OMP_THREAD_LIMIT --env OMP_WAIT_POLICY --compress\""
  #bash[15,1] <- "sort -u $SLURM_NODEFILE > uniqueNodes"
  #bash[16,1] <- "parallel --jobs 8 --sshloginfile uniqueNodes --workdir $PWD < REVEALERcommands.txt"
  #bash[17,1] <- ")"
  
  #write.table( bash,paste0(LCfolder,paste0(prefix,".bash", sep = ""), sep = ""), sep="\t",  col.names=FALSE,row.names = FALSE,na="",quote = FALSE)
  
}

RunREVEALER <- function(Target,
                        MicrobeList = F,
                        target.match = "positive",
                        Feature, 
                        seed.names = NULL,
                        exclude.features = NULL,
                        max.n.iter = 5,
                        pdf.output.file,
                        count.thres.low = NULL,
                        count.thres.high = NULL,
                        n.markers = 30,
                        locs.table.file = NULL){
  #args = commandArgs(trailingOnly=TRUE)
  source("REVEALER_library.R")
  
  if(MicrobeList == F){
    MicList <- names(BiMic)[3:length(BiMic[1,])]
  }
  else{
    MicList <- read.csv(MicrobeList,header=F,sep="\t")
    MicList <- MicList["V1"][,1]
  }
  
  for( i in MicList){
    REVEALER.v1(
      ds1                     = Target,  # Target profile dataset
      target.name             = i,       # Target profile
      target.match            = "negative",                    # Direction of match
      ds2                     = Feature,  # Features dataset
      seed.names              = seed.names,               # Seed(s) name(s)
      max.n.iter              = max.n.iter,                             # Maximum number of iterations
      pdf.output.file          = paste0(paste0("Revealer_REVEALER_",i),".pdf"),
      count.thres.low         = count.thres.low,                             # Filter out features with less than count.thres.low 1's
      count.thres.high        = count.thres.high,                            # Filter out features with more than count.thres.low 1's
      n.markers               = n.markers,                            # Number of top hits to show in heatmap for every iteration
    )
    break
  }
    
}
  
inputMaker<- function(gisticFile, mutFile,LCfolder){
  message(">>> making input for REVEALER, don't hurry! This will take a while so take a break : )")
  file<-read.delim(gisticFile,check.names = FALSE)
  fileMat<-format(file,trim=TRUE) #convert to matrix for faster processing
  fileMat<-fileMat[,-c(2,3)] # delete 2nd and 3rd columns
  mat<-matrix(,nrow=(2*nrow(fileMat)),ncol=ncol(fileMat)) #allocate empty matrix
  for (i in 1:nrow(fileMat)){
    # creates row labels for loci amplification or deletion
    paste(fileMat[i,1],"_AMP",sep="")->mat[i+(i-1),1]
    paste(fileMat[i,1],"_DEL",sep="")->mat[i*2,1]
    fileMat[i,1]<-0 #delete gene name
    as.matrix(which(fileMat[i,]>0))->indMatAmp #creates matrix of indices for samples with amplication
    mat[i+(i-1),indMatAmp]<-1
    as.matrix(which(fileMat[i,]<0))->indMatDel #for samples with copy number losses
    mat[2*i,indMatDel]<-1
  }
  
  #get names of column headers
  samples<-as.matrix(colnames(fileMat))
  
  #open mutation annotation file
  mutf<-read.delim(mutFile)
  mutMat<-as.matrix(mutf)
  mutMat<-mutMat[,c(1,16)] #only keep columns 1 and 16--Gene name and patient ID
  mat2<-matrix(,nrow=1,ncol=ncol(mat)) #initializes row for first gene
  mat2[1,1]<-paste(mutMat[1,1],"_MUT",sep="") #name first row with gene symbol
  n=1 #initializes counter for mat2 row index
  id<-substr(mutMat[1,2],1,12) #find patient ID from sample barcode
  ind<-grep(id,samples) #get column index of patient with mutation
  mat2[n,ind]<-1 #mark mutation is true for that patient
  for (i in 2:nrow(mutMat)){
    if (mutMat[i,1]== mutMat[i-1,1]){ #don't need to establish new line if same gene
      id<-substr(mutMat[i,2],1,12) #find patient ID from sample barcode
      ind<-grep(id,samples) #get column index of patient with mutation
      mat2[n,ind]<-1 #mark mutation is true for that patient
    } else {
      tempMat<-matrix(,nrow=1,ncol=ncol(mat)) # create one row for new gene
      tempMat[1,1]<-paste(mutMat[i,1],"_MUT",sep="") #label gene name
      id<-substr(mutMat[i,2],1,12) #find patient ID from sample barcode
      ind<-grep(id,samples) #get column of patient with mutation
      tempMat[1,ind]<-1 #mark mutation is true for that patient
      mat2<-rbind(mat2,tempMat) #add new row into exisiting matrix
      n=n+1 #move down one row
    }
  }
  
  #fill missing values with 0
  mat[is.na(mat)]<-0
  mat2[is.na(mat2)]<-0
  #combine CNV and mutation matrices
  mat3<-rbind(mat,mat2)
  #insert second column 
  newCol<-matrix("na",nrow=nrow(mat3),ncol=1) 
  rowName<-mat3[,1]
  mat3<-mat3[,-1]
  mat3<-cbind(rowName,newCol,mat3)
  #set up gct file format
  gct<-matrix(,nrow=3,ncol=ncol(mat3))
  gct[1,1]<-"#1.2"
  gct[2,1]<-nrow(mat3)
  gct[2,2]<-ncol(mat)-1
  gct[3,]<-c("NAME","Description",samples[-1])
  CNVdata<<-rbind(gct,mat3)
  
  write.table(CNVdata,paste0(LCfolder,"Features.txt"),sep="\t",quote=FALSE,na="",row.names=FALSE,col.names=FALSE)
  
}

ContamCorr <- function(df,out='./Contamination/',cutoff = 0.1,adjplot='nonadj') {
  dir.create(file.path(out), showWarnings = FALSE)
  dir.create(file.path(paste0(out),"PlateCorrection/"), showWarnings = FALSE)
  clin <- df$clinical
  headCat <- cli[names(df$vars)[1],]
  headCat <- rbind(headCat,cli["patient.samples.sample.portions.portion.analytes.analyte.aliquots.aliquot.plate_id",])
  rownames(headCat)[1]<-"Patients"
  CatTrans <- as_tibble(t(headCat), rownames = "row_names")
  CatTrans <-CatTrans[,2:ncol(CatTrans)]
  CatCombined <- merge(df$comb[1:length(row.names(df$BiMic)),],CatTrans, by="Patients", all=F,suffixes="Patients")
  CatCombined <- subset(CatCombined, select=-c(Type))
  CatCombined <- CatCombined %>% filter(Patients != "Patients")
  CatCombined <- CatCombined
  Kruskal(task = "PlateCorrection",out = paste0(out,"PlateCorrection/"),matinput = CatCombined,CatNum = ncol(CatTrans)-1,NumNum = ncol(df$BiMic)-2,cutoff = cutoff,adjplot=adjplot)
}