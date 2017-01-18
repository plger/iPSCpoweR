.checkPkg <- function(package){
  return(tryCatch(require(package, character.only=TRUE),error = function(e) FALSE))
}

#' initPar
#'
#' Initializes multithreading. This does not need to be called directly; the package's functions will call this.
#'
#' @param ncores The number of cores to use, default detectedCores()-1.
#'
#' @return A makeCluster object if parallel is supported and more than one core is used; otherwise FALSE.
#'
#' @export
initPar <- function(ncores=NULL){
    hasP <- .checkPkg("parallel")
    if(is.null(ncores) & hasP){
        library(parallel)
        ncores <- detectCores() - 1
    }else{
        if(ncores>1){
		if(!hasP){
			warning("Could not find the `parallel` package. Continuing with single thread.")
			ncores <- 1
		}else{
			library(parallel)
		}
	}
    }
    if(ncores>1)	return(makeCluster(ncores))
    return(FALSE)
}

#' getTxExpr
#'
#' Returns the transcript-level expression of the HipSci dataset.
#'
#' @param unit The expression unit in which to return the data; either "counts" (default), "fpkm", or "tpm".
#'
#' @export
getTxExpr <- function(unit="counts"){
    unit <- match.arg(unit, c("tpm","fpkm","counts"))
    data("hipsci_txCounts")
    if(unit != "counts") txCounts <- fpkm(txCounts)
    if(unit=="tpm") txCounts <- fpkm2tpm(txCounts)
    txCounts <- txCounts[which(rowSums(txCounts)>0),]
    return(txCounts)
}


#' getGeneExpr
#'
#' Returns the gene-level expression of the HipSci dataset.
#'
#' @param unit The expression unit in which to return the data; either "counts" (default), "fpkm", or "tpm".
#'
#' @export
getGeneExpr <- function(unit="counts", quiet=F){
    unit <- match.arg(unit, c("tpm","fpkm","counts"))
    txCounts <- getTxExpr(unit)
    data("tx2gene")
    if(!quiet)  message(paste("Aggregating transcript",unit,"to gene-level..."))
    gc <- aggregate(txCounts,by=list(gene=tx2gene[row.names(txCounts)]),FUN=sum)
    row.names(gc) <- gc$gene;
    gc$gene <- NULL
    return(gc)
}

#' Convert FPKM values to TPM values
#'
#' Convert FPKM values to transcripts per million (TPM).
#'
#' @param fpkm A dataframe, matrix or numeric vector of fpkm values.
#'
#' @return The corresponding TPM values.
#'
#' @export
fpkm2tpm <- function(fpkm){
  if(is.matrix(fpkm) | is.data.frame(fpkm)){
    for(i in 1:ncol(fpkm))    fpkm[,i] <- fpkm2tpm(as.numeric(fpkm[,i]))
    return(fpkm)
  }
  tpm <- fpkm/sum(fpkm)*10^6
  tpm[which(is.na(tpm))] <- 0
  return(tpm)
}

#' Convert counts to FPKM values
#'
#' Convert fragment counts to FPKM values, using the sum as library size. By default, effective length is used; set mean.frag.size to 0 to use total (given) length.
#'
#' @param counts A numeric dataframe or matrix, with transcripts as rows.
#' @param mean.frag.size A numeric value indicating the mean fragment size (default 220) used for effective length calculation.
#' @param removeShorter Logical; whether to remove transcripts whose length is shorter than the mean fragment size (default FALSE).
#'
#' @return A dataframe containing the corresponding FPKM values.
#'
#' @export
fpkm <- function(counts, mean.frag.size=500, removeShorter=F){
    data("txLengths")
    r <- intersect(names(txLengths),row.names(counts))
    if(length(r)==0)    stop("Could not find lengths for these transcripts. Is the data really transcript counts?")
    if(length(r) < nrow(counts))    message(paste("Could not find length for",nrow(counts)-length(r),"transcripts; these will be excluded."))
    txLengths <- txLengths[r]-mean.frag.size
    counts <- counts[r,]
    if(removeShorter) counts[which(txLengths < 0),] <- 0
    txLengths[which(txLengths < mean.frag.size)] <- mean.frag.size
    txLengths <- txLengths/1000
    for(i in 1:ncol(counts))    counts[,i] <- 10^6*counts[,i]/(txLengths*sum(counts[,i]))
    as.data.frame(counts)
}



#' aggByClone
#'
#' Aggregate expression data from technical replicates
#'
#' @param dat The expression data. If NULL (default), transcript-level count data are used.
#' @param removeFeeder Logical; whether to exclude samples from lines grown on feeder cells (default TRUE).
#' @param ag.fun Aggregation function; default `sum`
#'
#' @return A list with two slots: `dat` containing the aggregated data, and `annotation` containing the corresponding metadata.
#'
#' @export
aggByClone <- function(dat=NULL, removeFeeder=T, ag.fun=sum){
    data("hipsci_annotation")
    if(is.null(dat)){
        data("hipsci_txCounts")
        dat <- txCounts
    }else{
        annotation <- annotation[colnames(dat),]
    }
    if(removeFeeder){
        w <- which(!annotation$feeder)
        annotation <- annotation[w,]
        dat <- dat[,w]
    }
    ag <- aggregate(t(dat),by=list(clone=annotation$line),na.rm=T,FUN=ag.fun)
    row.names(ag) <- ag$clone
    ag$clone <- NULL;
    ag <- as.data.frame(t(ag))
    annotation <- annotation[!duplicated(annotation$line),]
    row.names(annotation) <- annotation$line
    return(list(annotation=annotation, dat=ag))
}


#' probSpurious
#'
#' Returns the probability of a given gene and/or foldchange in the permutation analyses.
#'
#' @param geneSymbol A HGNC gene symbol
#' @param log2FC The log2(foldchange) in your test (optional, but either geneSymbol or log2FC must be given)
#'
#' @return The probability of the gene and/or foldchange in the permutation analyses.
#'
#' @export
probSpurious <- function(geneSymbol=NULL, log2FC=NULL){
	if(is.null(geneSymbol) & is.null(log2FC)) stop("At least one of geneSymbol and log2FC must be given.")
	data("probSpurious")
	p1 <- ifelse( !is.null(geneSymbol) & geneSymbol %in% row.names(probGeneDE), probGeneDE[geneSymbol,"prob"], NA) 
	if(!is.na(p1))	print(paste("Gene's probability of differential expression across permutations:",format(p1, scientific=T)))
	p2 <- ifelse( is.null(log2FC), NA, probSpuriousFC(abs(log2FC)) )
	if(!is.na(p2))	print(paste("Probability of the foldchange in permutations:",format(p2, scientific=T)))
	if(is.na(p1) & is.na(p2))	message(paste("Could not find information for gene '",geneSymbol,"'",sep=""))
	return(p1,p2)
}

#' multiProbSpurious
#'
#' Annotates a data.frame of DEA results (such as produced by edgeR's `topTags` function) with the probability of each gene and foldchange in the DEA permutations.
#'
#' @param d A data.frame of differential expression results, but HGNC symbols as row.names and a column indicating the log2 foldchange.
#' @param FCfield The name of the column in `d` containing the log2(foldchange) (default `logFC`).
#' @param FDRfield Optional; the name of the column in `d` containing the FDR. If given, a column will be added with the maximum of the FDR and the probability of DE in the permuations. (default disabled)
#'
#' @return The data.frame d, with the addition columns probDEinPermutations and probFCinPermutations
#'
#' @export
multiProbSpurious <- function(d, FCfield="logFC", FDRfield=NULL){
	data("probSpurious")
	d <- as.data.frame(d)
	if(!(FCfield %in% colnames(d)))	stop("Could not find the foldchange column.")
	d$probDEinPermutations <- NA
	i <- intersect(row.names(d),row.names(probGeneDE))
	d[i,"probDEinPermutations"] <- probGeneDE[i,"prob"]
	d$probFCinPermutations <- probSpuriousFC(abs(d[i,FCfield]))
	if(!is.null(FDRfield)){
		if(!(FDRfield %in% colnames(d)))	stop("Could not find the column specified by the `FDRfield` argument.")
		w <- which(d$probFCinPermutations > 0.05 & d$probDEinPermutations > 0.01)
		d$maxFDR <- d[[FDRfield]]
		d[w,"maxFDR"] <- apply(cbind(d[w,FDRfield],d[w,"probDEinPermutations"]),1,na.rm=T,FUN=max)
	}
	return(d)
}


#' getSamplesInfo
#'
#' Returns an annotation of the samples used.
#'
#' @param removeFeeder Logical; whether to remove the samples cultivated on feeders, which should be exluded from the permutation analysis (Default TRUE).
#' @param doPlot Whether to plot summary information on the samples (default TRUE).
#'
#' @return A data.frame.
#'
#' @export
getSamplesInfo <- function(removeFeeder=T, doPlot=T){
	data("hipsci_annotation")
	if(removeFeeder){
		a <- annotation[which(!annotation$feeder),]
	}else{
		a <- annotation
	}
	a$Ethnicity <- as.factor(as.character(a$Ethnicity))
	if(!doPlot)	return(a[,c(1:4,6,8)])
	layout(matrix(1:4,nrow=2))
	nbClone <- table(aggregate(a$clone,by=list(indiv=a$individual),FUN=function(x){ length(unique(x))})$x)
	sex <- table(aggregate(a$sex,by=list(indiv=a$individual),FUN=function(x){ x[[1]] })$x)
	eth <- table(aggregate(a$Ethnicity,by=list(indiv=a$individual),FUN=function(x){ x[[1]] })$x)
	barplot(sex,ylab="Number of individuals",main="Sex")
	barplot(eth,ylab="Number of individuals",main="'Ethnicity'")
	barplot(nbClone,ylab="Number of individuals",main="iPSC clone per individual")
	return(a[,c(1:4,6,8)])
}